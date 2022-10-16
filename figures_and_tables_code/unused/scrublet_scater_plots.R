source("plotting.R")
library(dplyr)

PATH <- "/Volumes/scqc/code/scrublet/"

makeScatterPlots <- function(cells, clusters) {
  best_fit_line <- stat_smooth(method = "lm", col = "#C42126", se = FALSE, size = 1)
  
  for (sig in c("n_genes/doublet_score", "log2_n_genes/doublet_score")) {
    for (smp in unique(cells$Channel)) {
      if (sig == "n_genes/doublet_score") {
        metric.name <- "n_genes"
        sig <- "doublet_score"
      } else if (sig == "log2_n_genes/doublet_score") {
        metric.name <- "log2_n_genes"
        sig <- "doublet_score"
      }
      
      data <- data.frame(score=cells[[sig]], clusters=cells$louvain_labels, metric=cells[[metric.name]], demux_type=cells$demux_type, sample=cells$Channel) 
      data <- subset(data, sample == smp)
      colnames(data) <- c("score", "clusters", "metric", "demux_type", "sample") #rename data columns
      
      pearson <- round(cor(data$score, data$metric, method = "pearson"), 3)
      spearman <- round(cor(data$score, data$metric, method = "spearman"), 3)
      ttl <- ggtitle(paste0(task.name, "Sample: ", smp, "Pearson: ", pearson, "; Spearman: ", spearman))
      
      plot <- ggplot(data, aes(x = metric, y = score, color = demux_type)) + geom_point() + best_fit_line + labs(y=sig, x=metric.name) + ttl
      ggsave1(paste0(results.dir, metric.name, "-", sig, "_sample_", smp, ".pdf"),  plot, type = "u", n.clusters = 20)
    }
  }
}

tissues <- c("PanglaoDB/Liver",
             "PanglaoDB/Testis",
             "human_other/Olfactory_Epithelium",
             "human_other/krasnow_lung",
             "human_other/liver",
             "human_other/skin",
             "human_other/skin_cellbender",
             "tabula_muris/Heart_and_Aorta",
             "tabula_muris/Lung",
             "tabula_muris_smartseq2/Colon",
             "tabula_senis_30m_Large_Intestine")

for (tissue in tissues) {
  print(tissue)
  task.name <<- gsub("/", "-", tissue) 
  tissue <<- gsub("/", "_", tissue)
  
  results.dir <<- paste0(PATH, "scatter_plots/")
  dir.create(results.dir, showWarnings = FALSE)
  results.dir <<- paste0(results.dir, tissue, "/")
  dir.create(results.dir, showWarnings = FALSE)
  
  cells <- read.csv(paste0(PATH, tissue, "-cells.csv" ))
  cells$louvain_labels <- as.factor(cells$louvain_labels)
  cells[["log2_n_genes"]] <- log2(cells$n_genes)
  clusters <- read.csv(paste0(PATH, tissue, "-clusters.csv"))
  makeScatterPlots(cells, clusters)
}