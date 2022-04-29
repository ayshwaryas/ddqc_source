import pandas as pd
from ddqc_pipeline.config.config import OUTPUT_DIR_CLUSTER
import os


tissue = "Heart_and_Aorta"
path = OUTPUT_DIR_CLUSTER + f"task2/{tissue}/"

result_files = dict()
for folder in os.listdir(path):
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
print("\n\n\n\n\n")

method1 = "100genes-1.4-mad-2"
method2 = "50genes-1.4-mad-2"
method1_cells = result_files[method1]
method2_cells = result_files[method2]
for cluster in sorted(pd.unique(method2_cells.cluster_labels)):
    cluster_cells = method2_cells[method2_cells.cluster_labels == cluster]
    unique_cells = cluster_cells.index.difference(method1_cells.index)
    unique_pct = len(unique_cells) / len(cluster_cells.index) * 100
    if unique_pct > 80:
        print(round(unique_pct, 3))