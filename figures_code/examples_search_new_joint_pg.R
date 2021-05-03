source("figures_code/config.R")


results.dir <- OUTPUT.DIR
fout_ddqc <- "~/Downloads/unique_clusters_ddqc.csv"
fout_cutoff <- "~/Downloads/unique_clusters_cutoff.csv"
min.pct <- 0.85

write("project,tissue,cluster,annotation,cell_type,n_cells,percent_unique,qc_criteria,genes_median,mito_median,ribo_median,markers", fout_ddqc)
write("project,tissue,cluster,annotation,cell_type,n_cells,percent_unique,qc_criteria,genes_median,mito_median,ribo_median,markers", fout_cutoff)


for (project in dir(results.dir, recursive=FALSE)) {
  for (tissue in dir(paste0(results.dir, project, "/"), recursive=FALSE)) {
    print(paste(project, tissue, sep="_"))
    source.dir <- paste0(results.dir, project, "/", tissue, "/1.4-joint_clustering_old/")
    
    if (!file.exists(paste0(source.dir, "!cells.csv")) || !file.exists(paste0(source.dir, "!clusters.csv"))) {
      print(paste0(project, "-", tissue))
      next
    }
    
    if (tissue == "Heart_Circulation") {
      next
    }
    
    cells <- read.csv(paste0(source.dir, "!cells.csv"))
    rownames(cells) <- cells$barcodekey
    clusters <- read.csv(paste0(source.dir, "!clusters.csv"))
    rownames(clusters) <- clusters$cluster
    
    if (!("color" %in% colnames(cells))) {
      print(paste0(project, "-", tissue))
      next
    }
    
    for (cl in rownames(clusters)) {
      cluster.cells <- subset(cells, louvain_labels == cl)
      mad.only.cells <- rownames(subset(cluster.cells, color %in% c("MAD2 only")))
      cutoff.only.cells <- rownames(subset(cluster.cells, color %in% c("Cutoff only")))
      all.cells <- rownames(cluster.cells)
      
      if (length(all.cells) >= 30 && length(mad.only.cells) / length(all.cells) >= min.pct) {
        mrk <- strsplit(as.character(clusters[cl,]$markers), split=";")[[1]]
        qc.criteria.line <- ""
        if (clusters[cl,]$genes_median <= 200) {
          qc.criteria.line <- "n_genes"
        }
        if (clusters[cl,]$mito_median >= 10) {
          qc.criteria.line <- "percent_mito"
        }
        if (clusters[cl,]$genes_median <= 200 && clusters[cl,]$mito_median >= 10) {
          qc.criteria.line <- "n_genes and percent_mito"
        }
        write(paste(project, tissue, cl, as.character(clusters[cl,]$annotation), as.character(clusters[cl,]$cell_type), length(cluster.cells$barcodekey), round(length(mad.only.cells) / length(all.cells), 4), qc.criteria.line, as.character(clusters[cl,]$genes_median), as.character(clusters[cl,]$mito_median), as.character(clusters[cl,]$ribo_median), paste(mrk[1:min(25, length(mrk))], collapse=";", sep=""),  sep=","), file = fout_ddqc, append = TRUE)
      }
      
      if (length(all.cells) >= 30 && length(cutoff.only.cells) / length(all.cells) >= min.pct) {
        mrk <- strsplit(as.character(clusters[cl,]$markers), split=";")[[1]]
        qc.criteria.line <- ""
        if (clusters[cl,]$genes_median <= 200) {
          qc.criteria.line <- "n_genes"
        }
        if (clusters[cl,]$mito_median >= 10) {
          qc.criteria.line <- "percent_mito"
        }
        if (clusters[cl,]$genes_median <= 200 && clusters[cl,]$mito_median >= 10) {
          qc.criteria.line <- "n_genes and percent_mito"
        }
        write(paste(project, tissue, cl, as.character(clusters[cl,]$annotation), as.character(clusters[cl,]$cell_type), length(cluster.cells$barcodekey), round(length(cutoff.only.cells) / length(all.cells), 4), qc.criteria.line, as.character(clusters[cl,]$genes_median), as.character(clusters[cl,]$mito_median), as.character(clusters[cl,]$ribo_median), paste(mrk[1:min(25, length(mrk))], collapse=";", sep=""),  sep=","), file = fout_cutoff, append = TRUE)
      }
    }
  }
}
