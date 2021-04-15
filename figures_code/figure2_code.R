source("figures_code/config.R")


tissues <- c("tabula_muris/Heart_and_Aorta")


for (tissue in tissues) {
  task.name <<- gsub("/", "-", tissue)
  source.dir <<- paste0(OUTPUT.DIR, tissue, "/1.4-none-0/")
  results.dir <<- paste0(PATH, "figure2_plots/", task.name, "/")
  
  
  dir.create(paste0(PATH, "figure2_plots/"), showWarnings = FALSE)
  dir.create(paste0(results.dir), showWarnings = FALSE)
  
  file.copy(paste0(source.dir, "!cells.csv"), paste0(results.dir, "!cells.csv"))
  file.copy(paste0(source.dir, "!clusters.csv"), paste0(results.dir, "!clusters.csv"))
  file.copy(paste0(source.dir, "!markers.csv"), paste0(results.dir, "!markers.csv"))

  
  cells <- read.csv(paste0(results.dir, "!cells.csv"))
  cells$louvain_labels <- as.factor(cells$louvain_labels)
  clusters <- read.csv(paste0(results.dir, "!clusters.csv"))
  
  generatePlots(cells, clusters$cell_type, no.ct.labels=FALSE)
}