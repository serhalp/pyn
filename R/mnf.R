# vim: set filetype=r:

dyn.load (paste ("src/mnf_funcs", .Platform$dynlib.ext, sep = ""))

normalize.mnf <- function (object, verbose = TRUE, ...) {
    if (class (object) == "RGList") {
        nobject <- object

        if (verbose)
            cat ("Normalizing green channel...\n")
        nobject$G <- normalizeChannel (object$G, object$genes, verbose = verbose, ...)
        if (verbose)
            cat ("Normalizing red channel...\n")
        nobject$R <- normalizeChannel (object$R, object$genes, verbose = verbose, ...)
        if (verbose)
            cat ("All processing complete.\n")

        nobject
    } else {
        stop ("'object' must be of class RGList")
    }

    #object <- MA.RG (object)
    #asExprSet (object, idColumn = "PROBE_ID")
}

normalizeChannel <- function (channel, genes, ki = 2, kb = 20, summaryStatistic.i = "mean", summaryStatistic.b = "median", verbose = TRUE) {
    if (class (channel) == "matrix") {
        si <- rowify (summaryStatistic.i)
        sb <- rowify (summaryStatistic.b)
        probes <- genes$Status == "Probe"
        channel.probes <- channel[probes]

        if (verbose)
            cat ("\tComputing residuals...\n")
        res <- residuals.mnf (channel.probes, order (genes$PROBE_ID[probes]), ki, si)

        if (verbose)
            cat ("\tLocating grid neighbours...\n")
        neighbours <- knn.mnf (genes[probes, c ("X", "Y")], kb)

        if (verbose)
            cat ("\tCorrecting values...\n")
        nchannel <- channel
        ncells <- as.integer (prod (dim (neighbours)))
        mappedRes <- .C ("map_values", ncells, as.integer (neighbours), as.integer (res), as.integer (vector (length = ncells)), NAOK = TRUE, DUP = FALSE) [[4]]
        nchannel[probes] <- channel.probes - sb (matrix (mappedRes, nrow = length (res)))

        nchannel
    } else {
        stop ("'channel' must be a matrix")
    }
}

rowify <- function (fun) switch (fun, mean = rowMeans, median = rowMedians, min = rowMin, max = rowMax)

residuals.mnf <- function (channel, pos, k, sumStat) {
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

pcor <- function (v) { n <- length (v); cor (v[1:(n - 1)], v[2:n]) }
