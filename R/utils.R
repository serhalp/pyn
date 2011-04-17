# FIXME: Throw proper library imports in these functions.

downsample.matrix <- function (X, k) {
    group <- function (s) ((s - 1) * k + 1):(s * k)
    sapply (1:(nrow (X) / k), function (j)
        sapply (1:(ncol (X) / k), function (i)
            mean (X[group (i), group (j)], na.rm = TRUE)
        )
    )
}

exprs.ps <- function (e, indices, max.ps.size = 69)
    sapply (indices, function (ps) c (e[ps], rep (NA, max.ps.size - length (ps))))

f.test <- function (m) {
    library (limma)
    ft <- limma::bwss.matrix (m)
    return ((ft$bss / ft$bdf) / (ft$wss / ft$wdf))
}

res.as.matrix <- function (r, grid) {
    o <- as.double (rep (NA, max (grid[,1] + 1) * max (grid[,2] + 1)))
    null <- .C ("map_to_grid", as.integer (length (r)), as.double (r), as.integer (grid[,1]), as.integer (grid[,2]), o, DUP = FALSE, NAOK = TRUE)
    return (matrix (o, nrow = max (grid[,1]) + 1))
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

    plot (1:712, null, type = "p", cex = 0.3, ylim = range (c (m, mn), na.rm = TRUE), col = "red")
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

compute.degs <- function (batch, design, contrasts, which.contrasts = 1, n = Inf,
    sort.by = "logFC")
{
    library (limma)

    eset <- rma (batch, verbose = FALSE)
    fit <- lmFit (eset, design)
    cfit <- contrasts.fit (fit, contrasts)
    cfit <- eBayes (cfit)

    degs <- lapply (which.contrasts, function (coef)
        topTable (cfit, coef = coef, sort.by = sort.by, number = n)[, c ("ID", "logFC")])

    return (degs)
}

common.degs <- function (a, b, n = 100) {
    length (intersect (a$ID[1:n], b$ID[1:n])) / n
}

plot.common.degs <- function (a, b, n = 1:1000) {
    plot (n, sapply (n, common.degs, a = a[[1]], b = b[[1]]),
        type = "l", lty = 1, ylim = c (0, 1))

    if (length (a) > 1) {
        for (i in 2:length (a)) {
            lines (n, sapply (n, common.degs, a = a[[i]], b = b[[i]]),
                type = "l", lty = (i %% 6) + 1, ylim = c (0, 1))
        }
    }
}
