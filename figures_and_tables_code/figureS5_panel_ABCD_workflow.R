suppressPackageStartupMessages({
     library(SingleCellExperiment)
       library(scRNAseq)
       library(scater)
       library(flexmix)
       library(splines)
       library(BiocParallel)
       library(miQC)
   })

sce=readRDS("~/tm_heart_filter.rds")
sce=readRDS("~/Olfactory_Epithelium_filter.rds")

### Follow code for each input dataset
mt_genes <- grepl("^MT-",  rownames(sce))
feature_ctrls <- list(mito = rownames(sce)[mt_genes])
 feature_ctrls
sce <- addPerCellQC(sce, subsets = feature_ctrls,BPPARAM = MulticoreParam())
head(colData(sce))
metrics <- as.data.frame(colData(sce))
p <- ggplot(metrics,aes(x=detected,y=subsets_mito_percent)) + geom_point()
model <- mixtureModel(sce)
parameters(model)
head(posterior(model))
plotModel(sce, model)
plotFiltering(sce, model)

sce <- filterCells(sce, model)

### Save results
saveRDS(sce,"~/Olfactory_Epithelium_filter_miqc.rds")
saveRDS(sce,"~/tm_heart_filter_miqc.rds")

olf_miqc=readRDS("~/Olfactory_Epithelium_filter_miqc.rds")
tmh_miqc=readRDS("~/tm_heart_filter_miqc.rds")
