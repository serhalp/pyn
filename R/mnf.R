# vim:set filetype=r:

require (affy)
dyn.load (paste ("src/mnf_funcs", .Platform$dynlib.ext, sep = ""))

mnf <- function (batch, samples, interest = "probeset", bias = "grid", features.i = NULL, features.b = NULL, ...) {
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
    exprs.probes <- normalize.mnf (batch, features.i, features.b, ...)
    exprs (batch) <- exprs.probes

    # Compute a summarized gene expression value for each probeset, in each array
    exprs.genes <- summarize.mnf (batch)

    # Computes only the p-value (much faster than 't.test') for a two-sample
    # t-test with unequal sample sizes and equal variances
    t.test.p.value <- function (a, b) {
        na <- length (a)
        nb <- length (b)
        df <- na + nb - 2
        tstat <- ((mean (a) - mean (b)) / (sqrt (((na - 1) * var (a) + (nb - 1) * var (b)) / df) * sqrt (1 / na + 1 / nb)))
        return (2 * pt (-abs (tstat), df))
    }

    # Returns a data frame with log FC and t-test p-value for each gene
    num.genes <- length (featureNames (batch))
    diff.expr.stats <- function (pair) {
        arrays.first <- which (samples == pair[1])
        arrays.second <- which (samples == pair[2])

        # Compute p-value of two-sample t-test, for each gene
        p.values <- sapply (seq (num.genes),
            function (g) t.test.p.value (exprs.genes[g,arrays.first], exprs.genes[g,arrays.second]))

        return (data.frame (
            fold.change = rowMeans (exprs.genes[,arrays.first]) - rowMeans (exprs.genes[,arrays.second]),
            p.value = p.values,
            var = apply (exprs.genes[,arrays.first], 1, var) + apply (exprs.genes[,arrays.second], 1, var),
            row.names = featureNames (batch)))
    }

    # Compute log FC and t-test p-value for each gene, for each pair of samples
    num.samples <- length (unique (samples))
    # TODO: return variance
    return (lapply (combn (num.samples, 2, simplify = FALSE), diff.expr.stats))
}

summarize.mnf <- function (batch, summaryStatistic = "mean") {
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

pcor <- function (v, d = 1) { n <- length (v); cor (v[1:(n - d)], v[(1 + d):n]) }

pcors <- function (v, ds = 0:10)
    sapply (ds, pcor, v = v)

plotAutoCor <- function (v, ds = 0:10, ...) {
    plot (pcors (v, ds) ~ ds, type = "b", xlab = "Lag", ylab = "Genomic autocorrelation", ylim = c (0, 1), ...)
}

ssd.probeset.intra <- function (indices, values, sample.size = 200) {
    indices <- sample (indices, sample.size)
    ssd <- 0

    # FIXME: Rewrite without loops.
    for (ps in indices) {
        np <- length (ps)
        for (i in 1:(np - 1)) {
            p <- ps[i]
            o <- ps[(i + 1):np]
            ssd <- ssd + sum ((values[p] - values[o]) ^ 2)
        }
    }

    sqrt (ssd)
}

ssd.probeset.inter <- function (indices, values, sample.size = 200) {
    indices <- sample (indices, sample.size)
    ssd <- 0

    # FIXME: Rewrite without loops.
    for (i in 1:(sample.size - 1)) {
        ps <- indices[[i]]
        np <- length (ps)
        for (p in ps)
            for (o in indices[(i + 1):sample.size])
                ssd <- ssd + sum ((values[p] - values[o]) ^ 2)
    }

    sqrt (ssd)
}

ssd.probeset.ratio <- function (indices, values, sample.size.intra = 200, sample.size.inter = 200)
    ssd.probeset.intra (indices, values, sample.size.intra) / ssd.probeset.inter (indices, values, sample.size.inter)

vars.probeset <- function (batch) {
    indices <- indexProbes (batch, which = "pm")
    values <- exprs (batch)
    vars <- matrix (nrow = length (indices), ncol = ncol (values))

    for (i in 1:length (indices))
        vars[i,] <- apply (values[indices[[i]],], 2, var)

    vars
}

ftest.mnf <- function (batch, which = 1) {
    indices <- indexProbes (batch, which = "pm")
    values <- exprs (batch)[,which]
    global_mean <- mean (pm (batch)[,which]) # 'values' contains pm + mm + control

    # FIXME: This is despicable; rewrite without loops.
    inter_var <- 0
    intra_var <- 0
    for (ps in indices) {
        ps_mean <- mean (values[ps])
        inter_var <- inter_var + length (ps) * (ps_mean - global_mean) ^ 2
        intra_var <- intra_var + sum ((values[ps] - ps_mean) ^ 2)
    }

    return (((length (values) - length (indices)) / (length (indices) - 1)) * (inter_var / intra_var))
}
