from ddqc_pipeline.config.config import MITO_PREFIXES, RIBO_PREFIXES, OUTPUT_DIR
from ddqc_pipeline.filtering import initial_qc
from ddqc_pipeline.reading import get_project_info, read_tissue
from ddqc_pipeline.utils import save_to_csv


project = "tabula_muris"
tissue = "Heart_and_Aorta"
is_human, annotations = get_project_info(project=project, tissue=tissue)
adata = read_tissue(project=project, tissue=tissue, annotations=annotations)

mito_prefix = MITO_PREFIXES["human"] if is_human else MITO_PREFIXES["mouse"]
ribo_prefix = RIBO_PREFIXES["human"] if is_human else RIBO_PREFIXES["mouse"]
adata = initial_qc(adata, n_genes=100, percent_mito=100, mito_prefix=mito_prefix, ribo_prefix=ribo_prefix)
save_to_csv(adata, f"{OUTPUT_DIR}{project}/{tissue}/1.4-mad-2/!cells_initial_100genes_only.csv")