# vim:set filetype=r:

library (affy)
# FIXME: Path here should be relative to package root.
dyn.load (paste ("~/workspace/rmnf/src/mnf_funcs", .Platform$dynlib.ext, sep = ""))

mnf <- function (batch, samples, interest = "probeset", bias = "grid",
    features.i = NULL, features.b = NULL, verbose = TRUE, ...)
{
    if (is.null (features.i)) {
        features.i <- switch (interest,
            probeset = NULL, # FIXME: This is not very expressive...
            stop ("Need one of 'interest' or 'features.i'")
        )
    }

    if (is.null (features.b)) {
        features.b <- switch (bias,
            grid = indices2xy (seq (nrow (exprs (batch))), abatch = batch),
            none = NULL,
            stop ("Need one of 'bias' or 'features.b'")
        )
    }

    # Compute a normalized intensity for each probe in each array
    batch <- normalize.mnf (batch, features.i, features.b, verbose = verbose, ...)

    # Compute a summarized gene expression value for each probeset, in each array
    exprs.genes <- summarize.mnf (batch, verbose)

    # Computes only the p-value (much faster than 't.test') for a two-sample
    # t-test with unequal sample sizes and equal variances
    t.test.p.value <- function (a, b) {
        na <- length (a)
        nb <- length (b)
        df <- na + nb - 2
        tstat <- ((mean (a) - mean (b)) /
            (sqrt (((na - 1) * var (a) + (nb - 1) * var (b)) / df) *
                sqrt (1 / na + 1 / nb)
            )
        )
        return (2 * pt (-abs (tstat), df))
    }

    # Returns a data frame with replicate variance, log FC and t-test p-value for each gene
    num.genes <- length (featureNames (batch))
    diff.expr.stats <- function (pair) {
        arrays.first <- which (samples == pair[1])
        arrays.second <- which (samples == pair[2])

        # Compute p-value of two-sample t-test, for each gene
        p.values <- sapply (seq (num.genes), function (g)
            t.test.p.value (exprs.genes[g,arrays.first], exprs.genes[g,arrays.second]))

        return (data.frame (
            fold.change = rowMeans (exprs.genes[,arrays.first]) -
                rowMeans (exprs.genes[,arrays.second]),
            p.value = p.values,
            var = apply (exprs.genes[,arrays.first], 1, var) +
                apply (exprs.genes[,arrays.second], 1, var),
            row.names = featureNames (batch)))
    }

    # Compute some stats for each pair of samples
    if (verbose)
        cat ("Computing differential expression...\n")
    num.samples <- length (unique (samples))
    return (lapply (combn (num.samples, 2, simplify = FALSE), diff.expr.stats))
}

summarize.mnf <- function (batch, summary.stat = mean, verbose = TRUE) {
    if (verbose)
        cat ("Summarizing probeset intensities...\n")
    exprs.probes <- exprs (batch)
    # TODO: (Maybe) return ExpressionSet/eSet here rather than just a matrix.
    return (t (sapply (indexProbes (batch, which = "pm"), function (ps)
        apply (exprs.probes[ps,], 2, summary.stat))))
}

normalize.mnf <- function (batch, features.i, features.b, res.pre = NULL,
    dolog = TRUE, doexp = FALSE, use.median = FALSE, verbose = TRUE, ...)
{
    e <- exprs (batch)
    indices <- indexProbes (batch, which = "pm")

    if (dolog)
        e <- log2 (e)

    normalize.mnf.array <- function (a) {
        if (verbose)
            cat ("Normalizing array ", a, "...\n", sep = "")

        if (is.null (features.i)) {
            if (is.null (res.pre))
                res <- residuals.mnf.probeset (e[,a], indices, use.median)
            else
                res <- res.pre[,a]
        } else {
            # TODO: Refactor handling of bias/interest/features, to avoid nonsense like this.
            stop ("'probeset' is currently the only valid value for 'interest'.")
        }

        # Make sure non-pm probes are not considered grid neighbours.
        features.b.pm <- features.b
        features.b.pm[is.na (res),] <- NA
        return (normalizeChannel (e[,a], features.i = features.i,
            features.b = features.b.pm, res = res, verbose = verbose, ...)
        )
    }

    if (!is.null (features.b)) {
        # Apparently there is no other way to preserve the dimnames: matrix()
        # flattens matrix input; as.matrix() ignores 'dimnames' arg.  Argh.
        dimnames.bak <- dimnames (e)
        e <- sapply (seq (ncol (e)), normalize.mnf.array)
        dimnames (e) <- dimnames.bak
    }

    if (doexp)
        e <- 2 ^ e

    exprs (batch) <- e
    return (batch)
}

residuals.mnf.probeset <- function (values, indices, use.median = FALSE) {
    return (.Call ("affy_residuals", indices, values, as.logical (use.median)))
}

residuals.mnf.replicate <- function (values, samples) {
    residuals.mnf.replicate.sample <- function (s) {
        e <- values[,which (samples == s)]
        cols <- e - rowMeans (e)
        return (cols)
    }

    return (do.call (cbind, lapply (unique (samples), residuals.mnf.replicate.sample)))
}

# TODO: Refactor normalize.mnf + normalizeChannel (into one function?)
normalizeChannel <- function (channel, features.i, features.b, ki = 2, kb = 20,
    summary.stat.i = mean, summary.stat.b = mean, res = NULL, verbose = TRUE)
{
    if (!is.vector (channel) && !(is.matrix (channel) && ncol (channel) == 1))
        stop ("'channel' must be a vector or a 1-column matrix")

    if (is.null (res)) {
        if (verbose)
            cat ("    Computing residuals...\n")
        res <- residuals.mnf (channel, features.i, as.integer (ki), summary.stat.i)
    }

    if (verbose)
        cat ("    Locating neighbours in bias space...\n")
    neighbours <- knn.mnf (features.b, as.integer (kb))

    if (verbose)
        cat ("    Correcting values...\n")
    ncells <- as.integer (prod (dim (neighbours)))
    res.mapped <- matrix (.C ("map_values", ncells, as.integer (neighbours),
        as.double (res), as.double (rep (NA, ncells)), NAOK = TRUE, DUP = FALSE
    ) [[4]], nrow = length (res))

    b <- !is.na (res)
    channel[b] <- channel[b] - apply (res.mapped[b, , drop = FALSE], 1, summary.stat.b)
    return (channel)
}

residuals.mnf <- function (channel, pos, k, summary.stat) {
    if (is.null (pos))
        return (channel)

    neighbours <- knn.mnf (pos, k)
    ncells <- as.integer (prod (dim (neighbours)))
    mappedVals <- .C ("map_values", ncells, as.integer (neighbours),
        as.double (channel), as.double (rep (NA, ncells)), NAOK = TRUE, DUP = FALSE
    ) [[4]]
    channel - apply (matrix (mappedVals, nrow = length (channel)), 1, summary.stat)
}

knn.mnf <- function (v, k, ...) {
    if (is.vector (v) || (is.matrix (v) && dim (v) [2] == 1))
        knn.mnf.1D (v, k, ...)
    else if ((is.matrix (v) || is.data.frame (v)) && dim (v) [2] == 2)
        knn.mnf.2D (v[,1], v[,2], k, ...)
    else
        stop ("'v' must be a vector or a 1- or 2-column matrix")
}

knn.mnf.1D <- function (x, k) {
    n <- length (x)
    matrix (.C ("array_neighbours", as.integer (n), as.integer (x), as.integer (k),
        as.integer (rep (NA, n * k)), NAOK = TRUE, DUP = FALSE
    ) [[4]], nrow = n, ncol = k, byrow = TRUE)
}

knn.mnf.2D <- function (x, y, k) {
    n <- length (x)
    matrix (.C ("grid_neighbours", as.integer (n), as.integer (x), as.integer (y),
        as.integer (k), as.integer (rep (NA, n * k)), NAOK = TRUE, DUP = FALSE
    ) [[5]], nrow = n, ncol = k, byrow = TRUE)
}

# TODO: 'image.mnf.repres' and 'image.mnf.psres' can probably be combined
image.mnf.repres <- function (batch, samples, which = 1:length (batch),
    transfo = log2, shuffle = FALSE, col = pseudoPalette (
        low = "blue", high = "red", mid = "white"
    ), ...)
{
    num.probes <- nrow (exprs (batch))
    num.samples <- length (unique (samples))

    if (is.function (transfo))
        pm (batch) <- transfo (pm (batch))
    batch.res <- batch
    exprs (batch.res) <- matrix (NA, ncol = ncol (exprs (batch)),
        nrow = nrow (exprs (batch))
    )

    pm (batch.res) <- residuals.mnf.replicate (pm (batch), samples)
    image (batch.res[,which], transfo = NULL, col = col, ...)
    return (summary (pm (batch.res)))
}

image.mnf.psres <- function (batch, which = 1:length (batch), transfo = log2, 
    shuffle = FALSE, cutoff = 0.5, col = pseudoPalette (
        low = "blue", high = "red", mid = "white"
    ), ...)
{
    library (affyPLM)

    if (is.function (transfo))
        pm (batch) <- transfo (pm (batch))
    batch.res <- batch
    indices <- indexProbes (batch, which = "pm")

    image.mnf.psres.array <- function (e) {
        col <- residuals.mnf.probeset (exprs (batch)[,e], indices)
        if (shuffle)
            col[!is.na (col)] <- sample (col[!is.na (col)])
        col[abs (col) < cutoff] <- NA
        return (col)
    }

    exprs (batch.res) <- sapply (1:length (batch), image.mnf.psres.array)
    image (batch.res[,which], transfo = NULL, col = col, ...)
    return (summary (pm (batch.res)))
}
