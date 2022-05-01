from cell_classification.mad_filtering import filter_cells
import pandas as pd

# path = "Z:\\revisions\\revisions2\\adaptive_mad_threshold\\"
path = "C:\\Users\\misha\\Downloads\\adaptive_mad_threshold\\"
tissue = "heart"
cells = pd.read_csv(path + tissue + "_initial_clustering.csv", index_col=0)
output_dict = {"threshold": [], "cluster": [], "cells_filtered": [], "cells_filtered_pct": []}
for threshold in range(10, 36):
    threshold /= 10
    cells_filtered = filter_cells(cells, classification="cluster_labels", method="mad", threshold=threshold)
    for cl in sorted(cells.cluster_labels.unique()):
        cluster_cells = cells_filtered[cells_filtered.cluster_labels == cl]
        output_dict["threshold"].append(threshold)
        output_dict["cluster"].append(cl)
        output_dict["cells_filtered"].append(cluster_cells.cluster_labels_passed_qc.value_counts().get(False, 0))
        output_dict["cells_filtered_pct"].append(cluster_cells.cluster_labels_passed_qc.value_counts().get(False, 0) /
                                                 len(cluster_cells.index))


df = pd.DataFrame.from_dict(output_dict)
with open(path + tissue + "_filtered_cells_count.csv", "w") as file:
    file.write(df.to_csv())
