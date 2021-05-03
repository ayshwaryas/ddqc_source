source("figures_code/config.R")


results.dir <- OUTPUT.DIR
fout_ng_lower_75q <- paste0(PATH, "trends_n_genes_lower_75q.csv")
fout_ng_lower_max <- paste0(PATH, "trends_n_genes_lower_max.csv")
fout_mito_upper_25q <- paste0(PATH, "trends_mito_upper_25q.csv")
fout_mito_upper_min <- paste0(PATH, "trends_mito_upper_min.csv")


n_genes_lower_bound <- 200
mito_upper_bound <- 10

HEADER <- "project,tissue,cluster,annotation,cell_type,n_cells,genes_mean,genes_median,mito_mean,mito_median,ribo_mean,ribo_median,markers"
write(HEADER, fout_ng_lower_75q)
write(HEADER, fout_ng_lower_max)
write(HEADER, fout_mito_upper_25q)
write(HEADER, fout_mito_upper_min)

for (project in dir(results.dir, recursive=FALSE)) {
  for (tissue in dir(paste0(results.dir, project, "/"), recursive=FALSE)) {
    source.dir <- paste0(results.dir, project, "/", tissue, "/1.4-mad-2/")
    
    if (!file.exists(paste0(source.dir, "!cells.csv")) || !file.exists(paste0(source.dir, "!clusters.csv"))) {
      print(paste0(project, "-", tissue))
      next
    }
    
    cells <- read.csv(paste0(source.dir, "!cells.csv"))
    rownames(cells) <- cells$barcodekey
    clusters <- read.csv(paste0(source.dir, "!clusters.csv"))
    rownames(clusters) <- clusters$cluster
    
    for (cl in rownames(clusters)) {
      cluster.cells <- subset(cells, louvain_labels == cl)
      
      if (length(cluster.cells$barcodekey) >= 0) {
        mrk <- strsplit(as.character(clusters[cl,]$markers), split=";")[[1]]
        info <- paste(project, tissue, cl, as.character(clusters[cl,]$annotation), as.character(clusters[cl,]$cell_type), length(cluster.cells$barcodekey), as.character(clusters[cl,]$genes_mean), as.character(clusters[cl,]$genes_median), as.character(clusters[cl,]$mito_mean), as.character(clusters[cl,]$mito_median), as.character(clusters[cl,]$ribo_mean), as.character(clusters[cl,]$ribo_median), paste(mrk[1:min(25, length(mrk))], collapse=";", sep=""),  sep=",")
        if (summary(cluster.cells$n_genes)[[4]] <= n_genes_lower_bound) {
          write(info, file=fout_ng_lower_75q, append = TRUE)
        }
        if (summary(cluster.cells$n_genes)[[5]] >= n_genes_lower_bound) {
          write(info, file=fout_ng_lower_max, append = TRUE)
        }
        if (summary(cluster.cells$percent_mito)[[2]] >= mito_upper_bound) {
          write(info, file=fout_mito_upper_25q, append = TRUE)
        }
        if (summary(cluster.cells$percent_mito)[[1]]>= mito_upper_bound) {
          write(info, file=fout_mito_upper_min, append = TRUE)
        }
      }
    }
  }
}
