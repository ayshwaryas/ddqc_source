source("figures_code/config.R")

PATH <- OUTPUT.DIR

makeScatterPlots <- function(cells, clusters) {
  high.mito <- NULL
  for (row in 1:nrow(clusters)) {
    cl <- clusters[row, "cluster"]
    cluster.cells <- subset(cells, louvain_labels == cl)
    q1 <- summary(cluster.cells$percent_mito)[[2]]
    if (q1 > 10) {
      high.mito <- c(high.mito, cl)
      print(cl)
    }
  }
  
  results.dir <<- paste0(results.dir, "signature_plots/")
  do.call(file.remove, list(list.files(results.dir, full.names = TRUE)))
  dir.create(results.dir, showWarnings = FALSE)
  
  best_fit_line <- stat_smooth(method = "lm", col = "#C42126", se = FALSE, size = 1)
  
  for (sig in c("apoptosis", "apoptosis1", "apoptosis2", "apoptosis3", "dissociation", "er_stress", "go_mito_resp_chain_mouse", "go_mito_resp_chain_human", "mito_genes", "mito_ribo", "ribo_genes",
                "percent_ribo", "n_genes", "log2_n_genes", "n_genes/ribo", "log2_n_genes/ribo")) {
    
    if (sig == "n_genes/ribo") {
      metric.name <- "n_genes"
      sig <- "percent_ribo"
    } else if (sig == "log2_n_genes/ribo") {
      metric.name <- "log2_n_genes"
      sig <- "percent_ribo"
    } else {
      if (!(sig %in% colnames(cells))) {
        next
      }
      metric.name <- "percent_mito"
    }
    data <- data.frame(score=cells[[sig]], clusters=cells$louvain_labels, metric=cells[[metric.name]]) 
    colnames(data) <- c("score", "clusters", "metric") #rename data columns
    
    pearson <- round(cor(data$score, data$metric, method = "pearson"), 3)
    spearman <- round(cor(data$score, data$metric, method = "spearman"), 3)
    ttl <- ggtitle(paste0(task.name, " Pearson: ", pearson, "; Spearman: ", spearman))
    
    plot <- ggplot(data, aes(x = metric, y = score, color = clusters %in% high.mito)) + geom_point() + best_fit_line + labs(y=sig, x=metric.name) + ttl
    ggsave1(paste0(results.dir, metric.name, "-", sig, ".pdf"),  plot, type = "u", n.clusters = 20)
  }
}

tissues <- c("human_other/Adipose",
             "human_other/Kidney2",
             "human_other/Liver",
             "human_other/Krasnow_Lung",
             "PanglaoDB/Bone_Marrow",
             "PanglaoDB/Mammary_Gland",
             "PanglaoDB/Pancreatic_Islets",
             "PanglaoDB/Substantia_Nigra",
             "PanglaoDB/Testis",
             "tabula_muris_smartseq2/Bone_Marrow",
             "tabula_muris_smartseq2/Cerebellum",
             "tabula_muris_smartseq2/Colon",
             "tabula_muris_smartseq2/Heart_and_Aorta",
             "tabula_muris/Bladder",
             "tabula_muris/Heart_and_Aorta",
             "tabula_muris/Lung",
             "tabula_muris/Mammary_Gland",
             "tabula_muris/Tongue",
             "tabula_muris/Trachea",
             "human_other/Heart_Circulation")

for (tissue in tissues) {
  for (method in c("1.4-none-0", "1.4-mad-2")) {
    print(tissue)
    task.name <<- gsub("/", "-", tissue)
    
    results.dir <<- paste0(PATH, tissue, "/", method, "/")
    cells <- read.csv(paste0(results.dir, "!cells.csv"))
    cells$louvain_labels <- as.factor(cells$louvain_labels)
    cells[["log2_n_genes"]] <- log2(cells$n_genes)
    clusters <- read.csv(paste0(results.dir, "!clusters.csv"))
    makeScatterPlots(cells, clusters)
  }
}