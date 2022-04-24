import pandas as pd
from config.config import OUTPUT_DIR_CLUSTER
import os


tissue = "Lung"
method = "mad-2"
path = OUTPUT_DIR_CLUSTER + f"task4/{tissue}/"

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
print(" ", *[t for t in header], sep=",")
for i in range(len(header)):
    method1 = header[i]
    print(method1, end=",")
    for j in range(len(header)):
        method2 = header[j]
        intersection_size = len(result_files[method1].index.intersection(result_files[method2].index))
        print(intersection_size, end="")
        if j != len(header) - 1:
            print(",", end="")
    print()

