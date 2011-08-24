downsample.matrix <- function (X, k) {
    group <- function (s) ((s - 1) * k + 1):(s * k)
    sapply (1:(nrow (X) / k),
            function (j) sapply (1:(ncol (X) / k),
                                 function (i) mean (X[group (i), group (j)],
                                                    na.rm = TRUE)))
}

exprs.ps <- function (e, indices, max.ps.size = 69)
    sapply (indices, function (ps) c (e[ps], rep (NA, max.ps.size - length (ps))))

f.test <- function (m) {
    library (limma)
    ft <- limma::bwss.matrix (m)
    return ((ft$bss / ft$bdf) / (ft$wss / ft$wdf))
}

# FIXME: Rename to make non residual specific.
res.as.matrix <- function (r, grid) {
    o <- as.double (rep (NA, max (grid[, 1] + 1) * max (grid[, 2] + 1)))
    null <- .C ("map_to_grid", as.integer (length (r)), as.double (r),
                as.integer (grid[, 1]), as.integer (grid[, 2]), o, DUP = FALSE,
                NAOK = TRUE)
    return (matrix (o, nrow = max (grid[, 1]) + 1))
}

plot.res.stats <- function (m, mn) {
    null <- apply (m, 1, min, na.rm = TRUE)
    null[is.infinite (null)] <- NA
    null2 <- apply (m, 1, max, na.rm = TRUE)
    null2[is.infinite (null2)] <- NA
    null3 <- apply (m, 1, sd, na.rm = TRUE)

    nulln <- apply (mn, 1, min, na.rm = TRUE)
    nulln[is.infinite (nulln)] <- NA
    null2n <- apply (mn, 1, max, na.rm = TRUE)
    null2n[is.infinite (null2n)] <- NA
    null3n <- apply (mn, 1, sd, na.rm = TRUE)

    plot (1:712, null, type = "p", cex = 0.3,
          ylim = range (c (m, mn), na.rm = TRUE), col = "red")
    lines (spline (1:712, null), col = "red")
    lines (spline (1:712, rowMeans (m, na.rm = TRUE)), col = "blue")
    lines (spline (1:712, null2), col = "green")
    lines (spline (1:712, null3), col = "black")
    lines (spline (1:712, nulln), lty = "dashed", col = "red")
    lines (spline (1:712, rowMeans (mn, na.rm = TRUE)), lty = "dashed", col = "blue")
    lines (spline (1:712, null2n), lty = "dashed", col = "green")
    lines (spline (1:712, null3n), lty = "dashed", col = "black")
}

which.pm <- function (batch) {
    # TODO
}

compute.degs <- function (batch, design, contrasts, which.contrasts = NULL, n = Inf,
                          sort.by = "logFC", ...) {
    library (limma)

    eset <- rma (batch, verbose = FALSE)
    fit <- lmFit (eset, design)
    cfit <- contrasts.fit (fit, contrasts)
    cfit <- eBayes (cfit)

    if (is.null (which.contrasts))
        which.contrasts <- 1:ncol (contrasts)

    degs <- lapply (which.contrasts,
                    function (coef) topTable (cfit, coef = coef,
                                              sort.by = sort.by, number = n, ...))

    return (degs)
}

common.degs <- function (a, b, n = 100) {
    length (intersect (a$ID[1:n], b$ID[1:n])) / n
}

plot.common.degs <- function (a, b, n = 1:1000,
                              xlab = "Number of differentially expressed genes",
                              ylab = "% overlap between original and corrected", ...) {
    plot (n, 100 * sapply (n, common.degs, a = a[[1]], b = b[[1]]), type = "l",
          lty = 1, ylim = c (0, 100), xlab = xlab, ylab = ylab, ...)

    if (length (a) > 1) {
        for (i in 2:length (a)) {
            lines (n, 100 * sapply (n, common.degs, a = a[[i]], b = b[[i]]),
                   type = "l", lty = (i %% 6) + 1, ylim = c (0, 100))
        }
    }
}

cor.diag <- function (batch, pos, d = 1, dx = NULL, dy = NULL, res = FALSE) {
    if (res)
        indices <- indexProbes (batch, which = "pm")
    else
        mm (batch) <- NA

    if (is.null (dx) & is.null (dy))
        dx <- dy <- d

    return (apply (log2 (exprs (batch)), 2, function (a) {
        if (res)
            a <- residuals.mnf.probeset (a, indices)
        pm.matrix <- res.as.matrix (a, pos)
        pm.matrix.shifted <- pm.matrix[-(1:dx), -(1:dy)]
        pm.matrix <- pm.matrix[-((nrow (pm.matrix) - dx + 1):(nrow (pm.matrix))),
                               -((ncol (pm.matrix) - dy + 1):(ncol (pm.matrix)))]

        return (cor (as.vector (pm.matrix), as.vector (pm.matrix.shifted),
                     use = "complete.obs"))
    }))
}

cor.window <- function (batch, pos, res = FALSE) {
    if (res)
        indices <- indexProbes (batch, which = "pm")
    else
        mm (batch) <- NA

    return (apply (log2 (exprs (batch)), 2, function (a) {
        if (res)
            a <- residuals.mnf.probeset (a, indices)
        pm.matrix <- res.as.matrix (a, pos)
        n <- nrow (pm.matrix)
        m <- ncol (pm.matrix)
        pm.matrix.neighbours <- list (
            pm.matrix[2:(n - 1), 1:(m - 2)],
            pm.matrix[1:(n - 2), 2:(m - 1)],
            pm.matrix[2:(n - 1), 3:m],
            pm.matrix[3:n, 2:(m - 1)]
        )
        pm.matrix.neighbours.avg <- matrix (rowMeans (sapply (pm.matrix.neighbours,
                                                              as.vector),
                                                      na.rm = TRUE), nrow = n - 2)

        return (cor (as.vector (pm.matrix[2:(n - 1), 2:(m - 1)]),
                     as.vector (pm.matrix.neighbours.avg), use = "complete.obs"))
    }))
}

apply.res.probeset <- function (batch) {
    indices <- indexProbes (batch, which = "pm")
    res <- sapply (1:length (batch),
                   function (a) residuals.mnf.probeset (log2 (exprs (batch[, a])), indices))
    return (res)
}

hist.res.probeset <- function (batch, which = 1:length (batch),
                               res = apply.res.probeset (batch[, which]), main = "",
                               xlab = "Probe residuals", ylab = "Density", ...)
{
    apply (res, 2,
           function (r) {
               hist (res, breaks = 1000, main = main, xlab = xlab, ylab = ylab,
                     prob = TRUE, ...)
               abline (v = 0, col = "red")
           })
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
               yaxt = "n", ...)
    }

    lapply (which, image.mnf.psres.array)

    if (draw.legend) {
        levels <- seq (min (breaks), max (breaks), length = length (cols))
        image (levels, 1, matrix (levels, ncol = 1), col = cols, yaxt = "n",
               xlab = "", ylab = "", ...)
    }

    return (NULL)
}
