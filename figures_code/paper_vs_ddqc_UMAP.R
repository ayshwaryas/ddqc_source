source("plotting.R")


tm.dir <- "/Volumes/scqc/output_pg/tabula_muris/Lung/1.4-mad-2/"

cells <- read.csv(paste0(tm.dir, "!cells.csv"))
rownames(cells) <- cells$barcodekey
cells$louvain_labels <- as.factor(cells$louvain_labels)
clusters <- read.csv(paste0(tm.dir, "!clusters.csv"))
rownames(clusters) <- clusters$cluster
cell.types <- clusters$cell_type

percent.exclusive <- read.csv("~/Downloads/tm_ddqc_vs_paper.csv")$percent.exclusive

lbls <- NULL #create labels in the following format: cluster #, Panglao Cell Type \n annotated Cell Type
for (i in 1:length(cell.types)) {
  lbls <- c(lbls, paste0(i, " ", cell.types[i]))
}
ttl <<- ggtitle("ddqc_Lung")
names(lbls) <- 1:(length(lbls))

cells$category <- "Clusters shared betweed ddqc and paper"
#cells[percent.exclusive[cells$louvain_labels] >= 0.95, ]$category <- "Clusters at least 95% unique to ddqc"
cells[percent.exclusive[cells$louvain_labels] == 1, ]$category <- "Clusters fully unique to ddqc"

plot.cols <- c("Clusters shared betweed ddqc and paper" = "#DB6400", "Clusters at least 95% unique to ddqc" = "#16697A","Clusters fully unique to ddqc" = "#FFA62B")
data <- data.frame(UMAP1 = cells$umap1, UMAP2 = cells$umap2, cluster = cells$louvain_labels, category=cells$category, annotation=cells$annotations)
plot <- ggplot(data, aes(x=UMAP1, y=UMAP2, color=category, fill=cluster)) + geom_point(size = 1) + theme_umap + ttl + scale_fill_discrete(labels = lbls) + scale_color_manual(values = plot.cols)
for (cl in levels(data$cluster)) { #add cluster labels
  cluster.data <- subset(data, cluster == cl)
  plot <- plot + annotate("text", x = mean(cluster.data$UMAP1), y = mean(cluster.data$UMAP2), label=cl, size = 7, fontface=2)
}

n.clusters <- length(levels(cells$louvain_labels))
pdf(paste0("~/Downloads/exclusive_clusters_100.pdf"),width=10 + 2.5 * ceiling(n.clusters / 13),height=10,useDingbats = FALSE)
print(plot + no_bkg)
dev.off()
