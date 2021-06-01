tm.dir <- "/Volumes/scqc/output_pg/tabula_muris/"
fout <- paste0("~/Downloads/tm_ddqc_exclusive.csv")
write("tissue,cluster,cell_type,reason,n.genes.median,n.genes.75th, percent.mito.median, percent.mito.25th,markers", fout)


for (tissue in dir(tm.dir, recursive=FALSE)) {
  cells <- read.csv(paste0(tm.dir, tissue, "/1.4-mad-2/!cells.csv"))
  rownames(cells) <- cells$barcodekey
  clusters <- read.csv(paste0(tm.dir, tissue, "/1.4-mad-2/!clusters.csv"))
  rownames(clusters) <- clusters$cluster
  
  percent.mito.exclusive.cnt <- 0
  percent.mito.exclusive.text <- ""
  n.genes.exclusive.cnt <- 0
  n.genes.exclusive.text <- ""
  
  for (cl in 1:nrow(clusters)) {
    mrk <- strsplit(as.character(clusters[cl,]$markers), split=";")[[1]]
    ct <- clusters$cell_type[cl] 
    if (ct != "") {
      cells.cl <- subset(cells, louvain_labels == cl)
      percent.mito.median <- round(summary(cells.cl$percent_mito)[[3]], 3)
      percent.mito.25th <- round(summary(cells.cl$percent_mito)[[2]], 3)
      n.genes.median <- round(summary(cells.cl$n_genes)[[3]], 3)
      n.genes.75th <- round(summary(cells.cl$n_genes)[[4]], 3)
      
      if (percent.mito.25th >= 10) {
        percent.mito.exclusive.cnt <- percent.mito.exclusive.cnt + 1
        percent.mito.exclusive.text <- paste0(percent.mito.exclusive.text, ct, ";")
        write(paste(tissue, cl, clusters$cell_type[cl], "percent_mito", n.genes.median, n.genes.75th, percent.mito.median, percent.mito.25th, paste(mrk[1:min(25, length(mrk))], collapse=";", sep=""), sep=","), file=fout, append=TRUE)
      } else if (n.genes.75th <= 200) {
        n.genes.exclusive.cnt <- n.genes.exclusive.cnt + 1
        n.genes.exclusive.text <- paste0(n.genes.exclusive.text, ct, ";")
        write(paste(tissue, cl, clusters$cell_type[cl], "n_genes", n.genes.median, n.genes.75th, percent.mito.median, percent.mito.25th, paste(mrk[1:min(25, length(mrk))], collapse=";", sep=""), sep=","), file=fout, append=TRUE)
      }
    }
  }
  if (percent.mito.exclusive.text != "") {
    percent.mito.exclusive.text <- paste0(percent.mito.exclusive.cnt, " (", percent.mito.exclusive.text, ")")
  } else {
    percent.mito.exclusive.text <- percent.mito.exclusive.cnt
  }
  if (n.genes.exclusive.text != "") {
    n.genes.exclusive.text <- paste0(n.genes.exclusive.cnt, " (", n.genes.exclusive.text, ")")
  } else {
    n.genes.exclusive.text <- n.genes.exclusive.cnt
  }
  #write(paste(tissue, percent.mito.exclusive.text, n.genes.exclusive.text, sep=","), file=fout, append=TRUE)
}