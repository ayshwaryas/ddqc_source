library(Seurat)
source("r_code/filters.R")

mtx <- Read10X("/Volumes/scqc/data/mouse/tabula_muris/Heart_and_Aorta/Heart_and_Aorta-10X_P7_4/")
data <- CreateSeuratObject(counts = mtx)


data <- initialQC(data, basic.n.genes=100, basic.percent.mt=80)
df.qc <- ddqc.metrics(data)

data <- filterData(data, df.qc)
