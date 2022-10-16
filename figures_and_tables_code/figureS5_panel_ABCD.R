library(dplyr)
source("plotting.R")

convert.cell.names <- function(cell.names) {
  return(sapply(strsplit(cell.names, split="-"), last))
}

wd <- "/Volumes/scqc/"
tmh_miqc <- readRDS(paste0(wd, "data/Olfactory_Epithelium_filter_miqc.rds"))
miqc.cellnames <- rownames(colData(tmh_miqc))

path <- paste0(wd, "output_copies/output_pg11-14-21/human_other/Olfactory_Epithelium/")
ddqc <- read.csv(paste0(path, "1.4-mad-2/!cells.csv"))
ddqc$barcodekey <- convert.cell.names(ddqc$barcodekey)
all.cells <- read.csv(paste0(path, "1.4-none-0/!cells.csv"))
all.cells$barcodekey <- convert.cell.names(all.cells$barcodekey)
all.cells.clusters <- read.csv(paste0(path, "1.4-none-0/!clusters.csv"))

cells.union <- union(miqc.cellnames, ddqc$barcodekey)
cells.intersection <- intersect(miqc.cellnames, ddqc$barcodekey)
cells.miqc.unique <- setdiff(miqc.cellnames, ddqc$barcodekey)
cells.ddqc.unique <- setdiff(ddqc$barcodekey, miqc.cellnames)

all.cells <- subset(all.cells, barcodekey %in% cells.union)
all.cells$color <- "Neither"
all.cells[all.cells$barcodekey %in% cells.miqc.unique,]$color <- "MiQC only"
all.cells[all.cells$barcodekey %in% cells.ddqc.unique,]$color <- "ddqc only"
all.cells[all.cells$barcodekey %in% cells.intersection,]$color <- "Both methods"


no.ct.labels <- FALSE
lbls <- NULL #create labels in the following format: cluster #, Panglao Cell Type \n annotated Cell Type
for (i in 1:length(all.cells.clusters$cell_type)) {
  if (no.ct.labels) {
    lbls <- c(lbls, i)
  } else {
    lbls <- c(lbls, paste0(i, " ", all.cells.clusters$cell_type[i]))
  }
}
ttl <<- ggtitle("DDQC vs. MiQC, Olfactory Epithelium")
names(lbls) <- 1:(length(lbls))


plot.cols <- c("MiQC only" = "#000000", "ddqc only" = "#16697A", "Both methods" = "#FFA62B") #to keep consistent plot colors
all.cells$louvain_labels <- as.factor(all.cells$louvain_labels)
data <- data.frame(UMAP1 = all.cells$umap1, UMAP2 = all.cells$umap2, cluster = all.cells$louvain_labels, color=all.cells$color, annotation=all.cells$annotations)

fltplot <- ggplot(data, aes(x=UMAP1, y=UMAP2, color=color, fill=cluster)) + geom_point(size = 1) + 
  theme_umap + guides(colour = guide_legend(override.aes = list(size=2))) + 
  scale_fill_discrete(labels = lbls) + scale_color_manual(values = plot.cols) + ttl

for (cl in levels(all.cells$louvain_labels)) { #add cluster labels
  cluster.data <- subset(data, cluster == cl)
  fltplot <- fltplot + annotate("text", x = mean(cluster.data$UMAP1), y = mean(cluster.data$UMAP2), label=cl, size = 7, fontface=2)
}
ggsave1("~/Documents/plot.pdf", fltplot, n.clusters = length(unique(all.cells$louvain_labels)), type = "u")
