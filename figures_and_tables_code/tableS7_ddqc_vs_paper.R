tm.dir <- "/Volumes/scqc/output_pg/tabula_muris/"
fout <- paste0("~/Downloads/tm_ddqc_vs_paper.csv")
write("tissue,cluster,cell_type,number of cells, number of cells exclusive to ddqc,percent exclusive,n.genes.median,percent.mito.median,markers", fout)
annotations <- read.csv("/Volumes/scqc/data/mouse/tabula_muris/annotations_droplet.csv")


for (tissue in c("Lung")) {#dir(tm.dir, recursive=FALSE)) {
  ann.tissue <- annotations[annotations$tissue == tissue, ]
  
  cells <- read.csv(paste0(tm.dir, tissue, "/1.4-mad-2/!cells.csv"))
  rownames(cells) <- cells$barcodekey
  clusters <- read.csv(paste0(tm.dir, tissue, "/1.4-mad-2/!clusters.csv"))
  rownames(clusters) <- clusters$cluster
  
  exclusive.cnt <- 0
  exclusive.text <- ""
  
  for (cl in 1:nrow(clusters)) {
    mrk <- strsplit(as.character(clusters[cl,]$markers), split=";")[[1]]
    ct <- clusters$cell_type[cl] 
    cells.cl <- subset(cells, louvain_labels == cl)
    percent.mito.median <- round(summary(cells.cl$percent_mito)[[3]], 3)
    percent.mito.25th <- round(summary(cells.cl$percent_mito)[[2]], 3)
    n.genes.median <- round(summary(cells.cl$n_genes)[[3]], 3)
    n.genes.75th <- round(summary(cells.cl$n_genes)[[4]], 3)
    
    n.ddqc.exclusive.cells <- length(setdiff(sapply(rownames(cells.cl), gsub, pattern="-", replacement="_"), paste(tissue, "_", ann.tissue$cell, sep="")))
    pct <- round(n.ddqc.exclusive.cells / length(cells.cl$barcodekey), 3)
    if (pct >= 0) {
      exclusive.cnt <- exclusive.cnt + 1
      exclusive.text <- paste0(exclusive.text, ct, ";")
      write(paste(tissue, cl, clusters$cell_type[cl], length(cells.cl$barcodekey), n.ddqc.exclusive.cells, pct, n.genes.median, percent.mito.median, paste(mrk[1:min(25, length(mrk))], collapse=";", sep=""), sep=","), file=fout, append=TRUE)
    }
  }
  if (exclusive.text != "") {
    exclusive.text <- paste0(exclusive.cnt, " (", exclusive.text, ")")
  } else {
    exclusive.text <- percent.mito.exclusive.cnt
  }
  #write(paste(tissue, nrow(clusters), exclusive.text,  sep=","), file=fout, append=TRUE)
}
