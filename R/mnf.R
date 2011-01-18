# vim:set filetype=r:

require (affy)
dyn.load (paste ("src/mnf_funcs", .Platform$dynlib.ext, sep = ""))

mnf <- function (batch, samples, interest = "probeset", bias = "grid", features.i = NULL, features.b = NULL, verbose = TRUE, ...) {
    if (is.null (features.i)) {
        features.i <- switch (interest,
            genome = order (batch$genes$PROBE_ID), # This is MP_-specific
            probeset = NULL,
            stop ("Need one of 'interest' or 'features.i'")
        )
    }

    if (is.null (features.b)) {
        features.b <- switch (bias,
            grid = indices2xy (seq (nrow (exprs (batch))), abatch = batch),
            length = order (length (batch$genes[, "sequence"]), sample (1:length (batch))),
            none = NULL,
            stop ("Need one of 'bias' or 'features.b'")
        )
    }

    # Compute a normalized intensity for each probe in each array
    exprs.probes <- normalize.mnf (batch, features.i, features.b, verbose = verbose, ...)
    exprs (batch) <- exprs.probes

    # Compute a summarized gene expression value for each probeset, in each array
    exprs.genes <- summarize.mnf (batch, verbose)

    # Computes only the p-value (much faster than 't.test') for a two-sample
    # t-test with unequal sample sizes and equal variances
    t.test.p.value <- function (a, b) {
        na <- length (a)
        nb <- length (b)
        df <- na + nb - 2
        tstat <- ((mean (a) - mean (b)) / (sqrt (((na - 1) * var (a) + (nb - 1) * var (b)) / df) * sqrt (1 / na + 1 / nb)))
        return (2 * pt (-abs (tstat), df))
    }

    # Returns a data frame with replicate variance, log FC and t-test p-value for each gene
    num.genes <- length (featureNames (batch))
    diff.expr.stats <- function (pair) {
        arrays.first <- which (samples == pair[1])
        arrays.second <- which (samples == pair[2])

        # Compute p-value of two-sample t-test, for each gene
        p.values <- sapply (seq (num.genes),
            function (g) t.test.p.value (exprs.genes[g,arrays.first], exprs.genes[g,arrays.second]))

        return (data.frame (
            var = apply (exprs.genes[,arrays.first], 1, var) + apply (exprs.genes[,arrays.second], 1, var),
            fold.change = rowMeans (exprs.genes[,arrays.first]) - rowMeans (exprs.genes[,arrays.second]),
            p.value = p.values,
            row.names = featureNames (batch))
        )
    }

    # Compute some stats for each pair of samples
    if (verbose)
        cat ("Computing differential expression...\n")
    num.samples <- length (unique (samples))
    return (lapply (combn (num.samples, 2, simplify = FALSE), diff.expr.stats))
}

summarize.mnf <- function (batch, summaryStatistic = "mean", verbose = TRUE) {
    if (verbose)
        cat ("Summarizing probeset intensities...\n")
    stat <- colify (summaryStatistic)
    exprs.probes <- exprs (batch)
    return (t (sapply (indexProbes (batch, which = "pm"),
        function (ps) stat (exprs.probes[ps,]))))
}

normalize.mnf <- function (batch, features.i, features.b, dolog = TRUE, verbose = TRUE, ...) {
    exprs <- exprs (batch)
    indices <- indexProbes (batch, which = "pm")

    if (dolog)
        exprs <- log2 (exprs)

    normalize.mnf.array <- function (a) {
        if (verbose)
            cat ("Normalizing array ", a, "...\n", sep = "")

        if (is.null (features.i))
            res <- .Call ("affy_residuals", indices, exprs[,a])
        else
            stop ("'probeset' is currently the only valid value for 'interest'.")
        # Make sure non-pm probes are not considered grid neighbours
        features.b.pm <- features.b
        features.b.pm[is.na (res),] <- NA
        pmn <- normalizeChannel (exprs[,a], features.i = features.i, features.b = features.b.pm, res = res, verbose = verbose, ...)
        exprs[!is.na (pmn),a] <- pmn[!is.na (pmn)]
    }

    if (!is.null (features.b))
        return (sapply (seq (ncol (exprs)), normalize.mnf.array))
    else
        return (exprs)
}

normalizeChannel <- function (channel, features.i, features.b, ki = 2, kb = 20, summaryStatistic.i = "mean", summaryStatistic.b = "median", res = NULL, verbose = TRUE) {
    if (!is.vector (channel) && !(is.matrix (channel) && ncol (channel) == 1))
        stop ("'channel' must be a vector or a 1-column matrix")
    # For some reason, rowMeans and variants do not work on 1-column matrices
    if (ki <= 1 && kb <= 1)
        stop ("ki and kb must both be integers greater than one")

    si <- rowify (summaryStatistic.i)
    sb <- rowify (summaryStatistic.b)

    if (is.null (res)) {
        if (verbose)
            cat ("    Computing residuals...\n")
        res <- residuals.mnf (channel, features.i, as.integer (ki), si)
    }

    if (verbose)
        cat ("    Locating neighbours in bias space...\n")
    neighbours <- knn.mnf (features.b, as.integer (kb))

    if (verbose)
        cat ("    Correcting values...\n")
    ncells <- as.integer (prod (dim (neighbours)))
    mappedRes <- .C ("map_values", ncells, as.integer (neighbours), as.integer (res), as.integer (rep (NA, ncells)), NAOK = TRUE, DUP = FALSE) [[4]]

    b <- !is.na (res)
    channel[b] <- channel[b] - sb (matrix (mappedRes, nrow = length (res))[b,])
    return (channel)
}

rowify <- function (fun) switch (fun, mean = rowMeans, median = rowMedians, min = rowMin, max = rowMax)
colify <- function (fun) switch (fun, mean = colMeans, median = colMedians, min = colMin, max = colMax)

residuals.mnf <- function (channel, pos, k, sumStat) {
    if (is.null (pos))
        return (channel)

    neighbours <- knn.mnf (pos, k)
    ncells <- as.integer (prod (dim (neighbours)))
    mappedVals <- .C ("map_values", ncells, as.integer (neighbours), as.integer (channel), as.integer (rep (NA, ncells)), NAOK = TRUE, DUP = FALSE) [[4]]
    channel - sumStat (matrix (mappedVals, nrow = length (channel)))
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
    matrix (.C ("array_neighbours", as.integer (n), as.integer (x), as.integer (k), as.integer (rep (NA, n * k)), NAOK = TRUE, DUP = FALSE) [[4]], nrow = n, ncol = k, byrow = TRUE)
}

knn.mnf.2D <- function (x, y, k) {
    n <- length (x)
    matrix (.C ("grid_neighbours", as.integer (n), as.integer (x), as.integer (y), as.integer (k), as.integer (rep (NA, n * k)), NAOK = TRUE, DUP = FALSE) [[5]], nrow = n, ncol = k, byrow = TRUE)
}
