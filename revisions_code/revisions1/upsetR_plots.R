library(UpSetR)

path <- "/Volumes/scqc/output_pg/task4/"
tissue <- "Lung/"
wd <- paste0(path, tissue)

list.input <- list()
for (directory in list.dirs(wd, recursive = FALSE, full.names = FALSE)) {
  list.input[[directory]] <- read.csv(paste0(wd, directory, "/!cells.csv"))$barcodekey
}

#pdf("~/Documents/research/regevlab/PRIMES/ddqc/revisions/task1/20220116_heart_clustercomparison_upsetr.pdf")
upset(fromList(list.input),nsets=length(list.input), order.by = "freq", text.scale = 2)
#dev.off()
