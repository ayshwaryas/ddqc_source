import os

import pandas as pd
import scanpy as sc
import celltypist
from celltypist import models

models.download_models(force_update=False)

organism = "human"
project = "other"
tissue = "krasnow_lung"
directory = f"Z:\\data\\{organism}\\{project}\\{tissue}\\"

labels = None

for file in os.listdir(directory):
    if file.startswith("P"):
        adata = sc.read_10x_h5(f"Z:\\data\\{organism}\\{project}\\{tissue}\\{file}\\raw_feature_bc_matrix.h5")
        sc.pp.filter_cells(adata, min_genes=100)
        sc.pp.filter_genes(adata, min_cells=3)

        adata.var['mt'] = adata.var_names.str.startswith('MT-')
        sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
        adata = adata[adata.obs.pct_counts_mt < 80, :]

        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

        model = models.Model.load(model='Human_Lung_Atlas.pkl')
        predictions = celltypist.annotate(adata, model=model)
        current_labels = predictions.predicted_labels
        current_labels.index = [f"{file}-{t}" for t in current_labels.index]

        if labels is None:
            labels = predictions.predicted_labels
        else:
            labels = pd.concat([labels, predictions.predicted_labels])

with open("Z:\\revisions\\revisions2\\cell_classification\\krasnow_lung\\cell_typist.csv", "w") as file:
    file.write(labels.to_csv())
