import pandas as pd

path = "Z:\\revisions\\revisions2\\cell_classification\\"
tissue = "krasnow_lung"
final_ddqc_clustering_path = \
    "Z:\\output_copies\\output_pg_kuchroolab_copy\\human_other\\krasnow_lung\\1.4-mad-2\\!cells.csv"

# single_r = pd.read_csv(f"{path}{tissue}\\SingleR.csv", index_col=0)
# azimuth = pd.read_csv(f"{path}{tissue}\\azimuth.tsv", delimiter="\t", index_col=0)
cell_typist = pd.read_csv(f"{path}{tissue}\\cell_typist.csv", index_col=0)
classifications = pd.DataFrame(index=cell_typist.index, data={
    # "single_r": single_r["x"],
    # "azimuth": azimuth["predicted.celltype.l2"],
    "cell_typist": cell_typist["predicted_labels"],
})


ddqc_cell_data = pd.read_csv(f"{path}{tissue}\\!cells_initial.csv", index_col=0)
ddqc_cell_data.index = [t + "-1" for t in ddqc_cell_data.index]


summary = pd.DataFrame(index=classifications.index, data={
    "n_counts": ddqc_cell_data["n_counts"],
    "n_genes": ddqc_cell_data["n_genes"],
    "percent_mito": ddqc_cell_data["percent_mito"],
    "percent_ribo": ddqc_cell_data["percent_ribo"],

    "ddqc_cluster": ddqc_cell_data["cluster_labels"],
    # "single_r": classifications["single_r"],
    # "azimuth": classifications["azimuth"],
    "cell_typist": classifications["cell_typist"],
})

with open(f"{path}{tissue}\\summary.csv", "w") as file:
    file.write(summary.to_csv())
