source ("~/workspace/rmnf/R/mnf.R")
dump.mnf.degs <- function (batches, designs, contrasts, output.dir = ".", kb = 20) {
    gene.ranks <- function (fit) {
        gene.ranks.comp <- function (coef) {
            t <- topTable (fit, coef = coef, number = Inf, sort.by = "none")[,c ("ID", "logFC")]
            return (data.frame (ID = t$ID, rank = rank (t$logFC, ties.method = "random")))
        }

        l <- lapply (1:ncol (fit$contrasts), gene.ranks.comp)
        names (l) <- colnames (fit$contrasts)
        return (l)
    }

    get.degs.exp <- function (i, batch, design, contrast) {
        cat ("Batch ", i, "...\n")

        # Unnormalized pipeline
        cat ("\tk = 0\n")
        eset.rma <- rma (batch, verbose = FALSE)
        fit.rma <- lmFit (eset.rma, design)
        cfit.rma <- contrasts.fit (fit.rma, contrast)
        cfit.rma <- eBayes (cfit.rma)
        degs <- gene.ranks (cfit.rma)
        save (degs, file = paste (output.dir, paste (paste (i, "rma", sep = "_"), "rdata", sep = "."), sep = "/"))

        # Normalized pipeline with varying 'kb'
        g <- indices2xy (seq (nrow (exprs (batch))), abatch = batch)
        for (k in kb) {
            cat ("\tk =", k, "\n")
            batch.mnf <- normalize.mnf (batch, NULL, g, dolog = T, doexp = T, kb = k, verbose = FALSE)
            eset.mnf.rma <- rma (batch.mnf, verbose = FALSE)
            fit.mnf.rma <- lmFit (eset.mnf.rma, design)
            cfit.mnf.rma <- contrasts.fit (fit.mnf.rma, contrast)
            cfit.mnf.rma <- eBayes (cfit.mnf.rma)
            degs <- gene.ranks (cfit.mnf.rma)
            save (degs, file = paste (output.dir, paste (paste (i, "mnf-rma", k, sep = "_"), "rdata", sep = "."), sep = "/"))
        }
    }

    mapply (get.degs.exp, 1:length (batches), batches, designs, contrasts)
    return (NULL)
}
