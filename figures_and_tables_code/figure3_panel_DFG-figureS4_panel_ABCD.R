source("figures_code/config.R")


tissues <- c("tabula_muris/Heart_and_Aorta",
             "tabula_muris/Lung",
             "human_other/Krasnow_Lung",
             "human_other/Olfactory_Epithelium")

for (tissue in tissues) {
  task.name <<- gsub("/", "-", tissue)
  results.dir <<- paste0(OUTPUT.DIR, tissue, "/1.4-joint_clustering_old/")
  
  cells <- read.csv(paste0(results.dir, "!cells.csv"))
  cells$louvain_labels <- as.factor(cells$louvain_labels)
  clusters <- read.csv(paste0(results.dir, "!clusters.csv"))
  
  generatePlots(cells, clusters$cell_type,joint=TRUE)
}