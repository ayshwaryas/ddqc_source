from ddqc_pipeline.config.local_config import local


LOCAL = local  # needed for running outside the cluster

# PATHS
if not local:
    DATA_DIR = "/broad/kuchroolab/ayshwarya/scqc/data/"
    OUTPUT_DIR = "/broad/kuchroolab/ayshwarya/scqc/output_pg/"
else:  # for debug outside of cluster
    DATA_DIR = "Z:\\data\\"
    OUTPUT_DIR = ""
    OUTPUT_DIR_CLUSTER = "Z:\\output_pg\\"

# CLUSTERING
resolution = 1.4  # this resolution gives results closest to seurat

# FILTERING
basic_genes_filter = 100  # basic nGenes filter (performed for all methods)
basic_mito_filter = 80
MITO_PREFIXES = {"human": "MT-", "mouse": "mt-"}  # prefixes of mitochondrial genes
RIBO_PREFIXES = {"human": "^Rp[sl]\d", "mouse": "^Rp[sl]\d"}  # prefixes of ribosomal genes

# if true, filtering will be done for the selected metric
do_counts = True
do_genes = True
do_mito = True
do_ribo = False

# METHOD COMPARISON SCRIPTS
# none - no additional filtering; cutoff - min 200 genes, max 10% mito; outlier and mad - data driven methods
MC_METHODS = (("none", 0), ("cutoff", 10), ("mad", 2))
MC_TASKS_PER_TISSUE = len(MC_METHODS)
