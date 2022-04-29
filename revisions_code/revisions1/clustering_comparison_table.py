import os

import pandas
import pandas as pd

from ddqc_pipeline.config.config import OUTPUT_DIR_CLUSTER

tissue = "Lung"
method = "mad-2"
path = OUTPUT_DIR_CLUSTER + f"task1/{tissue}/"

total_cells = len(pandas.read_csv(path + "louvain-1.4-none-0/!cells.csv").barcodekey)
print(total_cells)

result_files = dict()
for folder in os.listdir(path):
    if method in folder:
        result_files[folder] = pd.read_csv(path + folder + "/!cells.csv", index_col="barcodekey")

common_cells = None
for cells in result_files.values():
    if common_cells is None:
        common_cells = cells.index
    else:
        common_cells = common_cells.intersection(cells.index)

all_except_leiden = None
for method, cells in result_files.items():
    if "leiden" not in method:
        if all_except_leiden is None:
            all_except_leiden = cells.index
        else:
            all_except_leiden = all_except_leiden.union(cells.index)

for method, cells in result_files.items():
    max_cluster_unique_pct = 0
    cluster_unique_number = 0
    for cluster in pd.unique(cells.cluster_labels):
        cluster_cells = cells[cells.cluster_labels == cluster]
        unique_cells = cluster_cells.index.difference(common_cells)
        unique_pct = len(unique_cells) / len(cluster_cells.index) * 100
        if max_cluster_unique_pct < unique_pct:
            max_cluster_unique_pct = unique_pct
            cluster_unique_number = cluster

    print(f"{method},{len(cells.index)},{round(max_cluster_unique_pct, 3)},{cluster_unique_number}")
print(f"Common: {len(common_cells)} cells")
print("\n\n\n")

header = sorted(list(result_files.keys()))
print(" ", *[t.split("-")[0] for t in header], sep=",")
for i in range(len(header)):
    method1 = header[i]
    print(method1.split("-")[0], end=",")
    for j in range(len(header)):
        method2 = header[j]
        intersection_size = len(result_files[method1].index.intersection(result_files[method2].index))
        print(intersection_size, end="")
        if j != len(header) - 1:
            print(",", end="")
    print()

leiden_cells = result_files["leiden-1.4-mad-2"]
louvain_cells = result_files["louvain-1.4-mad-2"]
leiden_unique = leiden_cells.loc[leiden_cells.index.difference(louvain_cells.index)]

louvain_filtered_counts = pandas.read_csv(path + "louvain-1.4-mad-2/!filtered_counts.csv")
louvain_filtered_genes = pandas.read_csv(path + "louvain-1.4-mad-2/!filtered_genes.csv")
louvain_filtered_mito = pandas.read_csv(path + "louvain-1.4-mad-2/!filtered_mito.csv")
leiden_unique["removal_reason"] = ["n_genes" if t in list(louvain_filtered_genes.barcodekey) else (
    "percent_ mito" if t in list(louvain_filtered_mito.barcodekey) else (
        "n_counts" if t in list(louvain_filtered_counts.barcodekey) else "ERROR")) for t in leiden_unique.index]

with open("/Users/michaelalperovich/Documents/tmp.csv", "w") as file:
    file.write(leiden_cells.loc[leiden_cells.index.difference(louvain_cells.index)].to_csv())
