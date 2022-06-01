library(ggplot2)
library(ggridges)
library(cowplot)
library(dplyr)

no_bkg <- theme(axis.line = element_line(colour = "black"),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(),
                panel.background = element_blank())
theme_horizontal <- theme(axis.text.x = element_text(angle = 45, size=15, hjust=1, face="bold"), 
                          axis.text.y = element_text(size=15), axis.title.y = element_text(size=15), 
                          legend.position="none", axis.title.x=element_blank(), plot.title = element_text(size = 20, face = "bold"))
theme_vertical <- theme(axis.text.x = element_text(size=15), axis.text.y = element_text(size=15, face="bold"), 
                        axis.title.x = element_text(size=15), legend.position="none", axis.title.y=element_blank(), 
                        plot.title = element_text(size = 20, face = "bold"))


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

generatePlotsByMetric <- function(obj, thresholds, lbls, metric.name.pg, metric.name, name.suffix) { #different plots for QC metrics
  labels_horizontal <- scale_x_discrete(labels=lbls)
  labels_vertical <- scale_y_discrete(labels=lbls)
  
  mean_plot <- stat_summary(fun="mean", geom="point", shape=23, fill="blue", size=3)
  
  if (metric.name.pg == "percent_mito") {
    axis_breaks_vertical <- scale_x_continuous(breaks=seq(0, 80, 5))
    vertical_line <- geom_vline(xintercept=10, color="red", size=0.5)
  } else {
    if (metric.name.pg == "n_genes") {
      vertical_line <- geom_vline(xintercept=log2(200), color="red", size=0.5)
    } else {
      vertical_line <- NULL
    }
    axis_breaks_vertical <- NULL
  }
  
  n.clusters <- length(levels(obj$cluster_labels))
  data <- data.frame(metric=obj[[metric.name.pg]], clusters=factor(obj$cluster_labels))
  colnames(data) <- c("metric", "clusters") #rename data column
  
  threshold1 <- thresholds[thresholds$threshold == 1, ]
  mad.barplot.data <- data.frame(mad_value=threshold1[[paste0(metric.name.pg, "_mad")]], clusters=factor(threshold1$cluster))
  colnames(mad.barplot.data) <- c("mad_value", "clusters") #rename data columns
  lines.data <- data.frame(MAD1=threshold1[[paste0(metric.name.pg, "_1MAD")]],
                           MAD2=threshold1[[paste0(metric.name.pg, "_2MAD")]],
                           MAD3=threshold1[[paste0(metric.name.pg, "_3MAD")]],
                           cluster=factor(threshold1$cluster))
  colnames(lines.data) <- c("MAD1", "MAD2", "MAD3", "cluster") #rename data columns
  
  if (metric.name.pg == "n_counts" || metric.name.pg == "n_genes") {
    # data$clusters = with(data, reorder(clusters, -log2(metric), mean)) #order data by cluster mean
    data$metric <- log2(data$metric)
    mad.barplot.data$mad_value <- log2(mad.barplot.data$mad_value)
    lines.data$MAD1 <- log2(lines.data$MAD1)
    lines.data$MAD2 <- log2(lines.data$MAD2)
    lines.data$MAD3 <- log2(lines.data$MAD3)
    
    axis_labels_vertical <- labs(x=paste0("log2(", metric.name, ")"))
    axis_labels_mad_horizontal <- labs(y=paste0("log2(MAD(", metric.name, "))"))
    axis_labels_mad_vertical <- labs(x=paste0("log2(MAD(", metric.name, "))"))
    name.suffix <- paste0(name.suffix, "_log.pdf")
    is.log2 = TRUE
  } else {
    # data$clusters = with(data, reorder(clusters, -metric, mean)) #order data by cluster mean
    axis_labels_vertical <- labs(x= metric.name)
    axis_labels_mad_horizontal <- labs(y=paste0("MAD(", metric.name, ")"))
    axis_labels_mad_vertical <- labs(x=paste0("MAD(", metric.name, ")"))
    name.suffix <- paste0(name.suffix, ".pdf")
    is.log2 = FALSE
  }
  
  joyplot <- ggplot(data, aes(x=metric, y=clusters)) + geom_density_ridges(
    jittered_points = TRUE,
    position = position_points_jitter(width = 0.05, height = 0),
    point_shape = '|', point_size = 3, point_alpha = 1, alpha = 0.7,
  ) + theme_vertical + mean_plot + labels_vertical + vertical_line + axis_breaks_vertical + ttl +
    geom_segment(data = lines.data, aes(x = MAD1, xend = MAD1, y = as.numeric(cluster), yend = as.numeric(cluster) + .9), color = "red") +
    geom_segment(data = lines.data, aes(x = MAD2, xend = MAD2, y = as.numeric(cluster), yend = as.numeric(cluster) + .9), color = "green") +
    geom_segment(data = lines.data, aes(x = MAD3, xend = MAD3, y = as.numeric(cluster), yend = as.numeric(cluster) + .9), color = "blue")
  
  mad.barplot <- ggplot(mad.barplot.data, aes(y=mad_value, x=clusters)) + geom_bar(stat="identity") + theme_horizontal + ttl
  mad.histogram <- ggplot(mad.barplot.data, aes(x=mad_value)) + geom_histogram() + theme_vertical + ttl
    
  ggsave1(filename=paste0(path, "density_", name.suffix), plot=joyplot + axis_labels_vertical, n.clusters = n.clusters, type = "v") 
  ggsave1(filename=paste0(path, "mad_barplot_", name.suffix), plot=mad.barplot + axis_labels_mad_horizontal, n.clusters = n.clusters, type = "h") 
  ggsave1(filename=paste0(path, "mad_histogram_", name.suffix), plot=mad.histogram + axis_labels_mad_vertical, n.clusters = n.clusters, type = "h") 
}

generateThresholdPlots <- function(thresholds) {
  mad_only <- thresholds[thresholds$threshold > 0,]
  outlier_only <- thresholds[thresholds$threshold == 0,]
  outlier_only$threshold = 1
  
  mad_sum <- mad_only %>% group_by(threshold) %>% summarise(cells_filtered = sum(cells_filtered))
  outlier_sum <- outlier_only %>% group_by(threshold) %>% summarise(cells_filtered = sum(cells_filtered))
  
  outlier_elevated <- outlier_only
  outlier_elevated$cells_filtered <- outlier_elevated$cells_filtered + 1
  outlier_elevated$cells_filtered_pct <- outlier_elevated$cells_filtered_pct + 1
  
  theshold_plot <- ggplot(mad_only, aes(threshold, cells_filtered)) + 
    geom_line() + facet_wrap(~cluster, ncol = 4, scales="free_y") + 
#    geom_hline(data = outlier_only, aes(yintercept = cells_filtered), color="red") + 
#    geom_text(data = outlier_elevated, label = outlier_only$cells_filtered, color="red") + 
    geom_vline(data = outlier_only, aes(xintercept = 2), color="blue")
  ggsave1(paste0(path, "threshold_vs_cells.pdf"), theshold_plot)
  
  theshold_pct_plot <- ggplot(mad_only, aes(threshold, cells_filtered_pct)) + geom_line() + 
    facet_wrap(~cluster, ncol = 4, scales="free_y") +
#    geom_hline(data = outlier_only, aes(yintercept = cells_filtered_pct), color="red") +
#    geom_text(data = outlier_elevated, label = outlier_only$cells_filtered_pct, color="red") +
    geom_vline(data = outlier_only, aes(xintercept = 2), color="blue")
  ggsave1(paste0(path, "threshold_vs_cells_pct.pdf"), theshold_pct_plot)
  
  theshold_aggregate_plot <- ggplot(mad_sum, aes(threshold, cells_filtered)) + geom_line() +
#    geom_hline(aes(yintercept = outlier_sum$cells_filtered), color="red") +
#    geom_text(data = outlier_sum, label = outlier_sum$cells_filtered, color="red") +
    geom_vline(data = outlier_sum, aes(xintercept = 2), color="blue")
  ggsave1(paste0(path, "threshold_vs_cells_sum.pdf"), theshold_aggregate_plot)
}


folder <- commandArgs(trailingOnly = TRUE)[1]
path <<- paste0("Z:\\revisions\\revisions2\\adaptive_mad_threshold\\", folder, "\\")
obj <- read.csv(paste0(path, "initial_clustering.csv"), row.names = 1)
# write.csv(obj, paste0(path, "initial_clustering.csv"))
thresholds <- read.csv(paste0(path, "filtered_cells_count.csv"), row.names = 1)
write.csv(thresholds, paste0(path, "filtered_cells_count.csv"))

path <<- paste0(path, "no_outlier_theshold_plots\\")

message("Making Plots")
lbls <- NULL #create labels in the following format: cluster #, Panglao Cell Type \n annotated Cell Type
for (i in 1:length(unique(obj$cluster_labels))) {
  lbls <- c(lbls, i)
}
ttl <<- ggtitle(folder)
names(lbls) <- 1:(length(lbls))

# generatePlotsByMetric(obj, thresholds, lbls, "n_counts", "Number of UMIS", "count")
# generatePlotsByMetric(obj, thresholds, lbls, "n_genes", "Number of Genes", "genes")
# generatePlotsByMetric(obj, thresholds, lbls, "percent_mito", "percent_mito", "mito")
# (obj, thresholds, lbls, "percent_ribo", "percent_ribo", "ribo")

generateThresholdPlots(thresholds)

