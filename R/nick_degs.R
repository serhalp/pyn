source ("/home/serhalp/workspace/rmnf/R/mnf.R")
dump.mnf.degs <- function (batches, designs, contrasts, output.dir = ".") {
    gene.ranks <- function (fit) {
        gene.ranks.comp <- function (coef) {
            t <- topTable (fit, coef = coef, number = Inf, sort.by = "none")[,c ("ID", "logFC")]
            return (data.frame (ID = t$ID, rank = order (t$logFC)))
        }

        l <- lapply (1:ncol (fit$contrasts), gene.ranks.comp)
        names (l) <- colnames (fit$contrasts)
        return (l)
    }

    get.degs.exp <- function (i, batch, design, contrast) {
        g <- indices2xy (seq (nrow (exprs (batch))), abatch = batch)
        batch.mnf <- normalize.mnf (batch, NULL, g, dolog = T, doexp = T, kb = 20)

        eset.rma <- rma (batch)
        eset.mnf.rma <- rma (batch.mnf)

        fit.rma <- lmFit (eset.rma, design)
        fit.mnf.rma <- lmFit (eset.mnf.rma, design)

        cfit.rma <- contrasts.fit (fit.rma, contrast)
        cfit.rma <- eBayes (cfit.rma)
        cfit.mnf.rma <- contrasts.fit (fit.mnf.rma, contrast)
        cfit.mnf.rma <- eBayes (cfit.mnf.rma)

        # Lousy R...
        degs.both <- lapply (list (cfit.rma, cfit.mnf.rma), gene.ranks)
        degs <- degs.both[[1]]
        save (degs, file = paste (output.dir, paste (paste (i, "rma", sep = "_"), "rdata", sep = "."), sep = "/"))
        degs <- degs.both[[2]]
        save (degs, file = paste (output.dir, paste (paste (i, "mnf-rma", sep = "_"), "rdata", sep = "."), sep = "/"))
    }

    mapply (get.degs.exp, 1:length (batches), batches, designs, contrasts)
    return (NULL)
}
