library(diptest)
library(LaplacesDemon)

tissues <- c("tabula_muris-Heart_and_Aorta",
           "tabula_muris-Lung",
           "human_other-Adipose",
           "human_other-Krasnow_Lung",
           "human_other-Olfactory_Epithelium",
           "human_other-Skin")


for (t in tissues) {
  project <- strsplit(t, "-")[[1]][1]
  tissue <- strsplit(t, "-")[[1]][2]
  cells.path <- paste0("Z:\\output_pg\\", project, "\\", tissue, "\\1.4-mad-2\\!cells_initial.csv")
  output.path <- paste0("Z:\\revisions\\revisions2\\adaptive_mad_threshold\\", project, "-", tissue, "\\cluster_modality.csv")
  cells <- read.csv(cells.path, row.names = 1)
  cells$cluster_labels <- as.factor(cells$cluster_labels)
  print(t)
  
  cluster <- NULL
  
  n_counts.diptest.p_val <- NULL
  n_counts.diptest.is_unimodal <- NULL
  n_counts.laplacesdemon.is_unimodal <- NULL
  
  n_genes.diptest.p_val <- NULL
  n_genes.diptest.is_unimodal <- NULL
  n_genes.laplacesdemon.is_unimodal <- NULL
  
  percent_mito.diptest.p_val <- NULL
  percent_mito.diptest.is_unimodal <- NULL
  percent_mito.laplacesdemon.is_unimodal <- NULL
  
  
  for (cl in levels(cells$cluster_labels)) {
    cells.cl <- cells[cells$cluster_labels == cl, ]
    cluster <- c(cluster, cl)
    
    n_counts.dt <- dip.test(cells.cl$n_counts)
    n_counts.diptest.p_val <- c(n_counts.diptest.p_val, n_counts.dt$p.value)
    n_counts.diptest.is_unimodal <- c(n_counts.diptest.is_unimodal, n_counts.dt$p.value > 0.05)
    n_counts.laplacesdemon.is_unimodal <- c(n_counts.laplacesdemon.is_unimodal, is.unimodal(cells.cl$n_counts))
    
    n_genes.dt <- dip.test(cells.cl$n_genes)
    n_genes.diptest.p_val <- c(n_genes.diptest.p_val, n_genes.dt$p.value)
    n_genes.diptest.is_unimodal <- c(n_genes.diptest.is_unimodal, n_genes.dt$p.value > 0.05)
    n_genes.laplacesdemon.is_unimodal <- c(n_genes.laplacesdemon.is_unimodal, is.unimodal(cells.cl$n_genes))
    
    percent_mito.dt <- dip.test(cells.cl$percent_mito)
    percent_mito.diptest.p_val <- c(percent_mito.diptest.p_val, percent_mito.dt$p.value)
    percent_mito.diptest.is_unimodal <- c(percent_mito.diptest.is_unimodal, percent_mito.dt$p.value > 0.05)
    percent_mito.laplacesdemon.is_unimodal <- c(percent_mito.laplacesdemon.is_unimodal, is.unimodal(cells.cl$percent_mito))
  }
  
  df <- data.frame(cluster=cluster,
    n_counts.diptest.p_val=n_counts.diptest.p_val, n_counts.diptest.is_unimodal=n_counts.diptest.is_unimodal, n_counts.laplacesdemon.is_unimodal=n_counts.laplacesdemon.is_unimodal,
    n_genes.diptest.p_val=n_genes.diptest.p_val, n_genes.diptest.is_unimodal=n_genes.diptest.is_unimodal, n_genes.laplacesdemon.is_unimodal=n_genes.laplacesdemon.is_unimodal,
    percent_mito.diptest.p_val=percent_mito.diptest.p_val, percent_mito.diptest.is_unimodal=percent_mito.diptest.is_unimodal, percent_mito.laplacesdemon.is_unimodal=percent_mito.laplacesdemon.is_unimodal
  )
  write.csv(df, output.path, row.names=FALSE)
}


