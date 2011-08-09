datasets.data <- function (celfile.path = "./data") {
    gse1400 <- ReadAffy (celfile.path = paste (celfile.path, "GSE1400/", sep = "/"))
    gse2189 <- ReadAffy (celfile.path = paste (celfile.path, "GSE2189/", sep = "/"))
    spikein <- ReadAffy (celfile.path = paste (celfile.path, "hgu133aspikein/", sep = "/"))
    batches <- list ("GSE1400" = gse1400, "GSE2189" = gse2189, "Spike-in" = spikein)
    return (batches)
}

datasets.data.corrected <- function (batches = datasets.data (), dolog = T, doexp = T, ...) {
    grid.133a <- indices2xy (seq (nrow (exprs (gse1400))), abatch = gse1400)
    grid.spikein <- indices2xy (seq (nrow (exprs (spikein))), abatch = spikein)
    grids <- list (grid.133a, grid.133a, grid.spikein)
    a <- list (5, 3, 16)
    batches.mnf <- mapply (function (batch, idx, g) normalize.mnf (batch[,idx], NULL, g,
        dolog = dolog, doexp = doexp, ...), batches, a, grids)
    return (batches.mnf)
}

datasets.contrasts <- function () {
    library (limma)
    samples <- list (sort (rep (1:2, 3)), sort (rep (1:6, 3)), c (sort (rep (10:14, 3)), sort (rep (1:9, 3))))
    designs <- lapply (samples, function (s) model.matrix (~ 0 + factor (s))
    colnames (designs[[1]]) <- c ("MEM", "CYT")
    colnames (designs[[2]]) <- c ("MGd_4h", "ctrl_4h", "MGd_12h", "ctrl_12h", "MGd_24h", "ctrl_24h")
    colnames (designs[[3]]) <- sapply (c (10:14, 1:9), function (i) paste ("exp", i, sep = "_"))
    contrasts <- list (makeContrasts (MEM - CYT, levels = designs[[1]]), makeContrasts (MGd_4h - ctrl_4h, MGd_12h - ctrl_12h, MGd_24h - ctrl_24h, levels = designs[[2]]), makeContrasts (exp_1 - exp_2, exp_2 - exp_3, exp_3 - exp_4, exp_11 - exp_12, exp_12 - exp_13, exp_13 - exp_14, levels = designs[[3]]))
}
