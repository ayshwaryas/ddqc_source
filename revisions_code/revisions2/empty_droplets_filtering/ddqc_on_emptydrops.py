import os

import numpy as np
import pandas as pd
import pegasus as pg

from ddqc_pipeline.config.config import DATA_DIR, MITO_PREFIXES, RIBO_PREFIXES, resolution, do_counts, do_genes, \
    do_mito, do_ribo
from ddqc_pipeline.filtering import initial_qc, filter_cells
from ddqc_pipeline.method_comparison import create_dirs
from ddqc_pipeline.reading import get_project_info, read_tissue


def initial_qc_emptydrops(adata, project, tissue, is_human, fdr_cutoff=0.01):
    species = "human" if is_human else "mouse"
    data_path = DATA_DIR + f"{species}/{project}/{tissue}/emptydrops/"

    emptydrops_results = None
    for sample in os.listdir(data_path):
        er = pd.read_csv(f"{data_path}/{sample}", sep="\t", index_col=0)
        er.index = [sample[:-4] + "-" + t[:-2] for t in er.index]
        if emptydrops_results is None:
            emptydrops_results = er.copy()
        else:
            emptydrops_results = pd.concat([emptydrops_results, er])

    adata["FDR"] = emptydrops_results.FDR
    adata = adata[adata.FDR <= fdr_cutoff]
    return adata


def main():
    project = "tabula_muris"
    tissue = "Heart_and_Aorta"
    method = "mad"
    param = 2

    is_human, annotations = get_project_info(project=project, tissue=tissue)
    adata = read_tissue(project, tissue, annotations)
    task_directory, task_name, results_dir = create_dirs(project, tissue, resolution, method, param)

    adata = initial_qc_emptydrops(adata, project, tissue, annotations)
    mito_prefix = MITO_PREFIXES["human"] if is_human else MITO_PREFIXES["mouse"]
    ribo_prefix = RIBO_PREFIXES["human"] if is_human else RIBO_PREFIXES["mouse"]
    adata = filter_cells(adata, resolution, method, param, basic_n_genes=0,
                         basic_percent_mito=100, mito_prefix=mito_prefix, ribo_prefix=ribo_prefix,
                         do_counts=do_counts, do_genes=do_genes, do_mito=do_mito, do_ribo=do_ribo,
                         record_path=results_dir, initial_result_name="initial_emptydrops")



