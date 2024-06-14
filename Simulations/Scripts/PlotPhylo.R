# Execute:
# Rscript --vanilla /home/torsten/Work/BDNN/ClusterScripts/PlotPhylo.R /home/torsten/Work/BDNN/BiSSE

library(ape)


getColorGradient <- function(X, Cols = c("red", "grey80", "blue")) {
  ColFunc <- colorRampPalette(Cols)
  M <- ceiling(100 * max(abs(range(X))))
  ColM <- ColFunc(2 * M + 1)
  ColIdx <- round(100 * X) + M
  return(ColM[ColIdx])
}


Args <- commandArgs(trailingOnly = TRUE)

Path <- Args
Replicate <- basename(Path)

Tree <- read.tree(file.path(Path, paste0(Replicate, '_phylo.tre')))

Traits <- read.table(file.path(Path, paste0(Replicate, '_traits_pvr.csv')),
                     header = TRUE, sep = '\t', row.names = 1)

TraitsOrdered <- Traits[match(Tree$tip.label, rownames(Traits)), ]

ColCat <- c('orange', 'dodgerblue')[TraitsOrdered$cat_trait_0 + 1]

ColCont <- getColorGradient(TraitsOrdered$cont_trait_0,
                              Cols = c("green4", "grey90", "magenta"))
CexCont <- TraitsOrdered$cont_trait_0
CexCont <- (CexCont - min(CexCont))
CexCont <- CexCont / max(CexCont)

CexPvr0 <- TraitsOrdered$pvr_0
CexPvr0 <- (CexPvr0 - min(CexPvr0))
CexPvr0 <- CexPvr0 / max(CexPvr0)
CexPvr1 <- TraitsOrdered$pvr_1
CexPvr1 <- (CexPvr1 - min(CexPvr1))
CexPvr1 <- CexPvr1 / max(CexPvr1)

pdf(file.path(Path, paste0(Replicate, '_phylo.pdf')),
    height = 6, width = 5, pointsize = 7)
par(mar = c(3, 0.5, 0.2, 0.2))
plot(Tree, show.tip.label = FALSE, root.edge = TRUE, x.lim = c(0, 36))
axisPhylo()
tiplabels(pch = 22, bg = ColCat, lwd = 0.5,
          frame = "n", cex = 0.8, offset = 0.05)
tiplabels(pch = 21, bg = 'green3', lwd = 0.5,
          frame = "n", cex = CexCont, offset = 0.5)
tiplabels(pch = 21, bg = 'grey', lwd = 0.5,
          frame = "n", cex = CexPvr0, offset = 1.0)
tiplabels(pch = 21, bg = 'grey', lwd = 0.5,
          frame = "n", cex = CexPvr1, offset = 1.5)
dev.off()

