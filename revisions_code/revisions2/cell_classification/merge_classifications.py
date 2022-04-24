import pandas as pd

path = "Z:\\cell_classification\\"
tissue = "pbmc"

single_r = pd.read_csv(path + tissue + "_SingleR.csv", index_col=0)
azimuth = pd.read_csv(path + tissue + "_azimuth.tsv", delimiter="\t", index_col=0)
cell_typist = pd.read_csv(path + tissue + "_cell_typist.csv", index_col=0)
classifications = pd.DataFrame(index=cell_typist.index, data={
    "single_r": single_r["x"],
    "azimuth": azimuth["predicted.celltype.l2"],
    "cell_typist": cell_typist["predicted_labels"],
})

metrics = pd.read_csv(path + tissue + "_metrics.csv")
metrics.index = [t + "-1" for t in metrics.barcodekey]
# ddqc_final_cluster_data = pd.read_csv("Z:\\output_pg\\task3\\misc\\pbmc_seurat\\1.4-mad-2\\!clusters.csv", index_col=1)
ddqc_final_cell_data = pd.read_csv("Z:\\output_pg\\task3\\misc\\pbmc_seurat\\1.4-mad-2\\!cells.csv")
retained_cells = [t[5:] + "-1" for t in ddqc_final_cell_data.barcodekey]
ddqc_cell_data = pd.read_csv(path + tissue + "_initial_clustering.csv")
ddqc_cell_data.index = [t + "-1" for t in ddqc_cell_data.barcodekey]
ddqc_cell_data["is_retained_by_ddqc"] = [t in retained_cells for t in ddqc_cell_data.index]
# ddqc_cell_data["cell_type"] = [ddqc_cluster_data.cell_type[t[1]["cluster_labels"]] for t in ddqc_cell_data.iterrows()]


summary = pd.DataFrame(index=classifications.index, data={
    "n_counts": metrics["n_counts"],
    "n_genes": metrics["n_genes"],
    "percent_mito": metrics["percent_mito"],
    "percent_ribo": metrics["percent_ribo"],

    "ddqc_cluster": ddqc_cell_data["cluster_labels"],
    "is_retained_by_ddqc": ddqc_cell_data["is_retained_by_ddqc"],
    # "ddqc_cell_type": ddqc_cell_data["cell_type"],

    "single_r": classifications["single_r"],
    "azimuth": classifications["azimuth"],
    "cell_typist": classifications["cell_typist"],
})

with open(path + tissue + "_summary.csv", "w") as file:
    file.write(summary.to_csv())
