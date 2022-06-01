library(ggplot2)
library(cowplot)

theme_umap <- theme(axis.text.x = element_text(size=15), axis.title.x = element_text(size=15),
                    axis.text.y = element_text(size=15), axis.title.y = element_text(size=15), 
                    plot.title = element_text(size = 20, face = "bold"), legend.title = element_text(size = 15), 
                    legend.text = element_text(size = 10))
theme_horizontal_with_legend <- theme(axis.text.x = element_text(angle = 45, size=15, hjust=1, face="bold"), 
                                      axis.text.y = element_text(size=15), axis.title.y = element_text(size=15), 
                                      axis.title.x=element_blank(), plot.title = element_text(size = 20, face = "bold"))
no_bkg <- theme(axis.line = element_line(colour = "black"),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(),
                panel.background = element_blank()) 


ggsave1 <- function(filename, plot, n.clusters=30, type="h") {  # custom ggsave function
  if (type == "h") {
    height = 10
    width = 14 / 30 * max(n.clusters, 30)
  }
  if (type == "v") {
    height = 14 / 30 * max(n.clusters, 30)
    width = 14
  }
  if (type == "u") {
    height = 10
    width = 10 + 2 * ceiling(n.clusters / 13)
  }
  
  ggsave(filename = filename, plot = plot + no_bkg, width = width, height = height) #saves plot with custom dimensions 
}

DimPlotCluster <- function(obj, lbls) { #DimPlot colored by cluster
  data <- data.frame(UMAP1 = obj$umap1, UMAP2 = obj$umap2, cluster = factor(obj$cluster_labels), other_cluster = factor(obj$other_cluster))
  plot <- ggplot(data, aes(x=UMAP1, y=UMAP2, color=other_cluster)) + geom_point(size = 1) + theme_umap + ttl + scale_fill_discrete(labels = lbls)
  for (cl in levels(obj$cluster_labels)) { #add cluster labels
    cluster.data <- subset(data, cluster == cl)
    plot <- plot + annotate("text", x = mean(cluster.data$UMAP1), y = mean(cluster.data$UMAP2), label=cl, size = 7, fontface=2)
  }
  return(plot)
}

makeBarplot <- function(obj, lbls) {
  data <- data.frame(UMAP1 = obj$umap1, UMAP2 = obj$umap2, cluster = obj$cluster_labels, other_cluster=obj$other_cluster)
  table.other_cluster <- NULL
  table.cluster <- NULL
  table.freq <- NULL
  for (cl in levels(data$cluster)) {
    data.cluster <- subset(data, cluster == cl)
    table.tmp <- as.data.frame(table(data.frame(other_cluster=factor(data.cluster$other_cluster), cluster = as.character(data.cluster$cluster))))
    table.tmp$Freq <- table.tmp$Freq / sum(table.tmp$Freq)
    table.other_cluster <- c(table.other_cluster, as.character(table.tmp$other_cluster))
    table.cluster <- c(table.cluster, as.integer(as.character(table.tmp$cluster)))
    table.freq <- c(table.freq, as.character(table.tmp$Freq))
  }
  data1 <- data.frame(other_cluster=table.other_cluster, cluster = as.factor(table.cluster), freq=as.double(as.character(table.freq)) * 100)
  freqplot <- ggplot(data1, aes(x=cluster, y=freq, fill=other_cluster)) + geom_bar(stat="identity") + theme_horizontal_with_legend + ttl + scale_x_discrete(labels=lbls)
  return(freqplot)
}


makePlots <- function(obj, results.dir, plot.title) {
  message("Making Plots")
  lbls <- NULL
  for (i in 1:length(unique(obj$cluster_labels))) {
    lbls <- c(lbls, i)
  }
  
  ttl <<- ggtitle(plot.title)
  names(lbls) <- 1:(length(lbls))
  n.clusters <- length(levels(obj$cluster_labels))
  ggsave1(filename=paste0(results.dir, plot.title, "_umap_clusters.pdf"), plot=DimPlotCluster(obj, lbls), n.clusters=n.clusters, type = "u")
  ggsave1(filename = paste0(results.dir, plot.title, "_barplot.pdf"), plot=makeBarplot(obj, lbls), n.clusters=n.clusters, type = "h")
}


folder <- "tabula_muris\\Lung"
path <<- paste0("Z:\\output_pg\\", folder, "\\1.4-mad-2\\")
cells.iqc <- read.csv(paste0(path, "!cells_initial.csv"), row.names = 1)
cells.emptydrops <- read.csv(paste0(path, "!cells_initial_emptydrops.csv"), row.names = 1)

cells.iqc$cluster_labels <- as.factor(cells.iqc$cluster_labels)
cells.iqc <- cells.iqc[cells.iqc$passed_qc == "True", ]
cells.iqc[["other_cluster"]] <- cells.emptydrops[rownames(cells.iqc), "cluster_labels"]

cells.emptydrops$cluster_labels <- as.factor(cells.emptydrops$cluster_labels)
cells.emptydrops <- cells.emptydrops[cells.emptydrops$passed_qc == "True", ]
cells.emptydrops[["other_cluster"]] <- cells.iqc[rownames(cells.emptydrops), "cluster_labels"]

makePlots(cells.iqc, path, plot.title = "iqc_colored_by_emptydrops")
makePlots(cells.emptydrops, path, plot.title = "emptydrops_colored_by_iqc")
