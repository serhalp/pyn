# vim: set filetype=r:

dyn.load (paste ("src/mnf_funcs", .Platform$dynlib.ext, sep = ""))

normalize <- function (object, verbose = TRUE, ...) {
    if (class (object) == "RGList") {
        nobject <- object

        if (verbose)
            cat ("Normalizing green channel...\n")
        nobject$R <- normalizeChannel (object$R, object$genes, verbose = verbose, ...)
        if (verbose)
            cat ("Normalizing red channel...\n")
        nobject$G <- normalizeChannel (object$G, object$genes, verbose = verbose, ...)
        if (verbose)
            cat ("All processing complete.\n")

        nobject
    } else {
        NA
    }

    #object <- MA.RG (object)
    #asExprSet (object, idColumn = "PROBE_ID")
}

normalizeChannel <- function (channel, genes, ki = 2, kb = 20, verbose = TRUE) {
    if (class (channel) == "matrix") {
        probes <- genes$Status == "Probe"
        channel.probes <- channel[probes]
        # TODO: sort by PROBE_ID

        if (verbose)
            cat ("\tComputing residuals...\n")
        res <- residuals.mnf (channel.probes, ki)

        if (verbose)
            cat ("\tLocating grid neighbours...\n")
        neighbours <- knn.mnf.2D (channel.probes, genes$X[probes], genes$Y[probes], kb)

        if (verbose)
            cat ("\tCorrecting values...\n")
        nchannel <- channel
        for (i in which (probes))
            nchannel[i] <- channel[i] - mean (res[neighbours[i,]])

        nchannel
    } else {
        NA
    }
}

residuals.mnf <- function (channel, k) {
    n <- length (channel)
    pred <- c (channel[2], rowMeans (cbind (channel[1:(n - 2)], channel[3:n])), channel[n - 1])
    channel - pred
}

knn.mnf.2D <- function (channel, x, y, k) {
    n <- length (channel)
    matrix (.C ("grid_neighbours", as.integer (n), as.integer (x), as.integer (y), as.integer (k), as.integer (vector (length = n * k)), NAOK = TRUE, DUP = FALSE) [[5]], nrow = n, ncol = k, byrow = TRUE)
}
