serhal.lemieux.2012 <- function (batch, ...) {
    return (normalize.pyn (batch, interest = "probeset", bias = "spatial", ...))
}

upton.lloyd.2005 <- function (batch, ...) {
    batch <- normalize.pyn (batch, res.pre = exprs (batch), bias = "spatial",
                            summary.stat.b = "min", kb = 8, do.log = FALSE,
                            do.exp = FALSE, ...)
    exprs (batch)[exprs (batch) < 0] <- 0
    return (batch)
}

normalize.pyn <- function (batch, interest = "probeset", bias = "spatial",
                           features.i = NULL, features.b = NULL, ki = 2, kb = 20,
                           summary.stat.i = "mean", summary.stat.b = "mean",
                           res.pre = NULL, do.log = TRUE, do.exp = TRUE,
                           verbose = TRUE) {
    # If no interest space provided, generate the requested one.
    if (is.null (features.i) && is.null (res.pre)) {
        features.i <- switch (interest,
                              probeset = NULL, # FIXME: This is not very expressive...
                              stop ("Need one of 'interest' or 'features.i'"))
    }

    # If no bias space provided, generate the requested one.
    if (is.null (features.b)) {
        features.b <- switch (bias,
                              spatial = indices2xy (seq (nrow (exprs (batch))), abatch = batch),
                              none = NULL,
                              stop ("Need one of 'bias' or 'features.b'"))
    }

    #stopifnot (summary.stat.b %in% c ("mean", "median", "min", "max"))

    e <- exprs (batch)
    indices <- indexProbes (batch, which = "pm")

    if (do.log)
        e <- log2 (e)

    normalize.pyn.array <- function (a) {
        if (verbose)
            cat ("Normalizing array ", a, " of ", length (batch), "...\n", sep = "")

        e.a <- e[, a]

        # TODO: Refactor handling of bias/interest/features, to avoid nonsense like this.
        if (is.null (res.pre)) {
            if (is.null (features.i)) {
                if (verbose)
                    cat ("    Computing residuals...\n")
                res <- residuals.pyn.probeset (e.a, indices, summary.stat.i)
            } else {
                stop ("'probeset' is currently the only valid value for 'interest'.")
            }
        } else {
            res <- res.pre[, a]
        }

        # Make sure non-pm probes are not considered grid neighbours.
        features.b.pm <- features.b
        features.b.pm[is.na (res), ] <- NA

        if (verbose)
            cat ("    Locating neighbours in bias space...\n")
        neighbours <- knn.pyn (features.b.pm, as.integer (kb))

        if (verbose)
            cat ("    Correcting values...\n")
        ncells <- as.integer (prod (dim (neighbours)))
        res.mapped <- matrix (.C ("map_values", ncells, as.integer (neighbours),
                                  as.double (res), as.double (rep (NA, ncells)),
                                  NAOK = TRUE, DUP = FALSE) [[4]],
                              nrow = length (res))

        b <- !is.na (res)
        e.a[b] <- e.a[b] - apply (res.mapped[b, , drop = FALSE], 1, summary.stat.b)
        return (e.a)
    }

    if (!is.null (features.b)) {
        # Apparently there is no other way to preserve the dimnames: matrix()
        # flattens matrix input; as.matrix() ignores 'dimnames' arg.  Argh.
        dimnames.bak <- dimnames (e)
        e <- sapply (seq (ncol (e)), normalize.pyn.array)
        dimnames (e) <- dimnames.bak
    } else {
        warning ("Warning: not actually normalizing arrays, as no bias space was specified.")
    }

    if (do.exp)
        e <- 2 ^ e

    exprs (batch) <- e
    return (batch)
}

pyn <- normalize.pyn

residuals.pyn.probeset <- function (values, indices, summary.stat = "mean") {
    stat.id <- switch (summary.stat,
                       mean = 0,
                       median = 1,
                       min = 2,
                       max = 3,
                       stop ("'summary.stat.i' must be 'mean', 'median', 'min' or 'max'"))

    return (.Call ("affy_residuals", indices, values, as.integer (stat.id)))
}

residuals.pyn.replicate <- function (values, samples) {
    residuals.pyn.replicate.sample <- function (s) {
        e <- values[, which (samples == s)]
        cols <- e - rowMeans (e)
        return (cols)
    }

    return (do.call (cbind, lapply (unique (samples), residuals.pyn.replicate.sample)))
}

knn.pyn <- function (v, k, ...) {
    if (is.vector (v) || (is.matrix (v) && dim (v) [2] == 1))
        knn.pyn.1D (v, k, ...)
    else if ((is.matrix (v) || is.data.frame (v)) && dim (v) [2] == 2)
        knn.pyn.2D (v[, 1], v[, 2], k, ...)
    else
        stop ("'v' must be a vector or a 1- or 2-column matrix")
}

knn.pyn.1D <- function (x, k) {
    n <- length (x)
    matrix (.C ("array_neighbours", as.integer (n), as.integer (x), as.integer (k),
        as.integer (rep (NA, n * k)), NAOK = TRUE, DUP = FALSE
    ) [[4]], nrow = n, ncol = k, byrow = TRUE)
}

knn.pyn.2D <- function (x, y, k) {
    n <- length (x)
    matrix (.C ("grid_neighbours", as.integer (n), as.integer (x), as.integer (y),
        as.integer (k), as.integer (rep (NA, n * k)), NAOK = TRUE, DUP = FALSE
    ) [[5]], nrow = n, ncol = k, byrow = TRUE)
}
