# vim:set filetype=r:

library (affy)
# FIXME: Path here should be relative to package root.
dyn.load (paste ("~/workspace/rmnf/src/mnf_funcs", .Platform$dynlib.ext, sep = ""))

# TODO: Consider default value of 'do.exp'.
normalize.mnf <- function (batch, interest = "probeset", bias = "spatial",
                           features.i = NULL, features.b = NULL, ki = 2, kb = 20,
                           summary.stat.i = mean, summary.stat.b = mean,
                           res.pre = NULL, do.log = TRUE, do.exp = FALSE,
                           verbose = TRUE, ...) {
    # If no interest space provided, generate the requested one.
    if (is.null (features.i)) {
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

    # Validate the summary statistics.
    use.median <- switch (summary.stat.i,
                          mean = FALSE,
                          median = TRUE,
                          stop ("'summary.stat.i' must be either 'mean' or 'median'"))

    e <- exprs (batch)
    indices <- indexProbes (batch, which = "pm")

    if (do.log)
        e <- log2 (e)

    normalize.mnf.array <- function (a) {
        if (verbose)
            cat ("Normalizing array ", a, " of ", length (batch), "...\n", sep = "")

        e.a <- e[, a]

        # TODO: Refactor handling of bias/interest/features, to avoid nonsense like this.
        if (is.null (res.pre)) {
            if (is.null (features.i)) {
                if (verbose)
                    cat ("    Computing residuals...\n")
                res <- residuals.mnf.probeset (e.a, indices, use.median)
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
        neighbours <- knn.mnf (features.b, as.integer (kb))

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
        e <- sapply (seq (ncol (e)), normalize.mnf.array)
        dimnames (e) <- dimnames.bak
    } else {
        warning ("Warning: not actually normalizing arrays, as no bias space was specified.")

    if (do.exp)
        e <- 2 ^ e

    exprs (batch) <- e
    return (batch)
}

residuals.mnf.probeset <- function (values, indices, use.median = FALSE) {
    return (.Call ("affy_residuals", indices, values, as.logical (use.median)))
}

residuals.mnf.replicate <- function (values, samples) {
    residuals.mnf.replicate.sample <- function (s) {
        e <- values[, which (samples == s)]
        cols <- e - rowMeans (e)
        return (cols)
    }

    return (do.call (cbind, lapply (unique (samples), residuals.mnf.replicate.sample)))
}

knn.mnf <- function (v, k, ...) {
    if (is.vector (v) || (is.matrix (v) && dim (v) [2] == 1))
        knn.mnf.1D (v, k, ...)
    else if ((is.matrix (v) || is.data.frame (v)) && dim (v) [2] == 2)
        knn.mnf.2D (v[, 1], v[, 2], k, ...)
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
                              transfo = log2, shuffle = FALSE,
                              col = rainbow (50), ...) {
    num.probes <- nrow (exprs (batch))
    num.samples <- length (unique (samples))
    if (shuffle)
        cat ("Apparently, shuffle is not actually implemented.")

    if (is.function (transfo))
        pm (batch) <- transfo (pm (batch))
    batch.res <- batch
    exprs (batch.res) <- matrix (NA, ncol = ncol (exprs (batch)),
                                 nrow = nrow (exprs (batch)))

    pm (batch.res) <- residuals.mnf.replicate (pm (batch), samples)
    image (batch.res[, which], transfo = NULL, col = col, ...)
    return (summary (pm (batch.res)))
}

image.mnf.psres <- function (batch, grid, which = 1:length (batch),
                             transfo = log2, draw.legend = TRUE,
                             shuffle = FALSE, ...) {
    breaks <- c (-10, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 10)
    cols <- c (rgb (0, 0, seq (1, 0.5, -0.1)), rgb (seq (0.5, 1, 0.1), 0, 0))

    if (is.function (transfo))
        pm (batch) <- transfo (pm (batch))
    batch.res <- batch
    indices <- indexProbes (batch, which = "pm")

    image.mnf.psres.array <- function (e) {
        col <- residuals.mnf.probeset (exprs (batch)[, e], indices)
        if (shuffle)
            col[!is.na (col)] <- sample (col[!is.na (col)])
        m <- res.as.matrix (col, grid)
        image (m[, !is.na (m[225, ])], col = cols, breaks = breaks, xaxt = "n",
               yaxt = "n")
    }

    lapply (which, image.mnf.psres.array)

    if (draw.legend) {
        levels <- seq (min (breaks), max (breaks), length = length (cols))
        image (levels, 1, matrix (levels, ncol = 1), col = cols, yaxt = "n",
               xlab = "", ylab = "")
    }

    return (NULL)
}
