# vim:set filetype=r:

require (affy)
dyn.load (paste ("src/mnf_funcs", .Platform$dynlib.ext, sep = ""))

normalize.mnf <- function (object, channel = "both", interest = "genome", bias = "grid", features.i = NULL, features.b = NULL, probe.filter = function (v) rep (TRUE, length (v)), verbose = TRUE, ...) {
    if (is.null (features.i)) {
        features.i <- switch (interest,
            genome = order (RG$genes$PROBE_ID), # This is MP_-specific
            none = NULL,
            stop ("Need one of 'interest' or 'features.i'")
        )
    }

    if (is.null (features.b)) {
        features.b <- switch (bias,
            grid = RG$genes[, c ("X", "Y")],
            length = order (length (RG$genes[, "sequence"]), sample (1:length (RG))),
            none = NULL,
            stop ("Need one of 'bias' or 'features.b'")
        )
    }

    return (switch (class (object),
        RGList = normalize.mnf.RGList (object, channel, features.i, features.b, probe.filter, verbose, ...),
        AffyBatch = normalize.mnf.AffyBatch (object, channel, features.i, features.b, probe.filter, verbose, ...),
        stop ("'object' must be of class RGList or AffyBatch")
    ))
}

normalize.mnf.RGList <- function (rg, channel, features.i, features.b, probe.filter, verbose, ...) {
    #probes <- rg$genes$Status == "Probe"
    probes <- probe.filter (rg$genes)

    for (e in 1:ncol (rg)) {
        if (verbose)
            cat ("Normalizing array ", e, "...\n", sep =  "")

        if (channel %in% c ("green", "both")) {
            if (verbose)
                cat ("  Green channel...\n")
            rg$G[probes] <- normalizeChannel (rg$G[probes,e], features.i = subset (features.i, probes), features.b = subset (features.b, probes), verbose = verbose, ...)
        }

        if (channel %in% c ("red", "both")) {
            if (verbose)
                cat ("  Red channel...\n")
            rg$R[probes] <- normalizeChannel (rg$R[probes,e], features.i = subset (features.i, probes), features.b = subset (features.b, probes), verbose = verbose, ...)
        }
    }

    if (verbose)
        cat ("All processing complete.\n")
    rg
    #rg <- MA.RG (rg)
    #asExprSet (rg, idColumn = "PROBE_ID")
}

normalize.mnf.AffyBatch <- function (batch, channel, features.i, features.b, probe.filter, verbose, ...) {
    npm <- matrix (nrow = nrow (pm (batch)), ncol = ncol (pm (batch)))
    nmm <- matrix (nrow = nrow (mm (batch)), ncol = ncol (mm (batch)))

    for (e in 1:length (batch)) {
        if (verbose)
            cat ("Normalizing array ", e, "...\n", sep =  "")

        if (channel %in% c ("pm", "both")) {
            if (verbose)
                cat ("  PM...\n")
            npm[,e] <- normalizeChannel (pm (batch[,e]), features.i = features.i, features.b = features.b, verbose = verbose, ...)
        }

        if (channel %in% c ("mm", "both")) {
            if (verbose)
                cat ("  MM...\n")
            nmm[,e] <- normalizeChannel (mm (batch[,e]), features.i = features.i, features.b = features.b, verbose = verbose, ...)
        }
    }

    if (channel %in% c ("pm", "both"))
        pm (batch) <- npm
    if (channel %in% c ("mm", "both"))
        mm (batch) <- nmm

    if (verbose)
        cat ("All processing complete.\n")
    batch
}

normalizeChannel <- function (channel, features.i, features.b, ki = 2, kb = 20, summaryStatistic.i = "mean", summaryStatistic.b = "median", verbose = TRUE) {
    if (is.vector (channel) || (is.matrix (channel) && ncol (channel) == 1)) {
        si <- rowify (summaryStatistic.i)
        sb <- rowify (summaryStatistic.b)

        if (verbose)
            cat ("    Computing residuals...\n")
        res <- residuals.mnf (channel, features.i, ki, si)

        if (verbose)
            cat ("    Locating neighbours in bias space...\n")
        neighbours <- knn.mnf (features.b, kb)

        if (verbose)
            cat ("    Correcting values...\n")
        ncells <- as.integer (prod (dim (neighbours)))
        mappedRes <- .C ("map_values", ncells, as.integer (neighbours), as.integer (res), as.integer (vector (length = ncells)), NAOK = TRUE, DUP = FALSE) [[4]]
        channel - sb (matrix (mappedRes, nrow = length (res)))
    } else {
        stop ("'channel' must be a vector or a 1-column matrix")
    }
}

rowify <- function (fun) switch (fun, mean = rowMeans, median = rowMedians, min = rowMin, max = rowMax)

residuals.mnf <- function (channel, pos, k, sumStat) {
    if (is.null (pos))
        return (channel)

    neighbours <- knn.mnf (pos, k)
    ncells <- as.integer (prod (dim (neighbours)))
    mappedVals <- .C ("map_values", ncells, as.integer (neighbours), as.integer (channel), as.integer (vector (length = ncells)), NAOK = TRUE, DUP = FALSE) [[4]]
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
    matrix (.C ("array_neighbours", as.integer (n), as.integer (x), as.integer (k), as.integer (vector (length = n * k)), NAOK = TRUE, DUP = FALSE) [[4]], nrow = n, ncol = k, byrow = TRUE)
}

knn.mnf.2D <- function (x, y, k) {
    n <- length (x)
    matrix (.C ("grid_neighbours", as.integer (n), as.integer (x), as.integer (y), as.integer (k), as.integer (vector (length = n * k)), NAOK = TRUE, DUP = FALSE) [[5]], nrow = n, ncol = k, byrow = TRUE)
}

pcor <- function (v, d = 1) { n <- length (v); cor (v[1:(n - d)], v[(1 + d):n]) }

pcors <- function (v, ds = 0:10) {
    cors <- vector (length = length (ds))
    for (i in 1:length (ds))
        cors[i] <- pcor (v, ds[i])
    cors
}

plotAutoCor <- function (v, ds = 0:10, ...) {
    plot (pcors (v, ds) ~ ds, type = "b", xlab = "Lag", ylab = "Genomic autocorrelation", ylim = c (0, 1), ...)
}

ssd.probeset.intra <- function (indices, values, sample.size = 200) {
    indices <- sample (indices, sample.size)
    ssd <- 0
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
