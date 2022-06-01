import pathlib

from revisions_code.revisions2.cell_classification.mad_filtering import filter_cells
import numpy as np
import pandas as pd

FILE_PATH = pathlib.Path(__file__).parent.resolve()

tissues = ["tabula_muris-Heart_and_Aorta",
           "tabula_muris-Lung",
           "human_other-Adipose",
           "human_other-Krasnow_Lung",
           "human_other-Olfactory_Epithelium",
           "human_other-Skin"]

for project, tissue in [t.split("-") for t in tissues]:
    cells_path = f"Z:\\output_pg\\{project}\\{tissue}\\1.4-mad-2\\!cells_initial.csv"
    output_path = pathlib.Path(f"Z:\\revisions\\revisions2\\adaptive_mad_threshold\\{project}-{tissue}\\")

    cells = pd.read_csv(cells_path + "", index_col=0)
    output_dict = {"threshold": [], "cluster": [], "cells_filtered": [], "cells_filtered_pct": []}
    metric_names = ["n_counts", "n_genes", "percent_mito"]
    for mn in metric_names:
        output_dict[mn + "_median"] = []
        output_dict[mn + "_mad"] = []
        output_dict[mn + "_sd"] = []
        output_dict[mn + "_mad_to_sd_ratio"] = []
        for i in range(1, 4):
            output_dict[f"{mn}_{i}MAD"] = []

    for threshold in range(10, 37):
        if threshold == 36:  # run outlier
            threshold = 0
            cells_filtered = filter_cells(cells, classification="cluster_labels", method="outlier", threshold=threshold)
        else:
            threshold /= 10
            cells_filtered = filter_cells(cells, classification="cluster_labels", method="mad", threshold=threshold,
                                          return_full=True)
        for cl in sorted(cells.cluster_labels.unique()):
            cluster_cells = cells_filtered[cells_filtered.cluster_labels == cl]
            n_filtered_cells = cluster_cells.cluster_labels_passed_qc.value_counts().get(False, 0)
            output_dict["threshold"].append(threshold)
            output_dict["cluster"].append(cl)
            output_dict["cells_filtered"].append(n_filtered_cells)
            output_dict["cells_filtered_pct"].append(round(n_filtered_cells / len(cluster_cells.index) * 100, 3))
            for mn in metric_names:
                if threshold == 0:
                    output_dict[mn + "_median"].append(0)
                    output_dict[mn + "_mad"].append(0)
                    output_dict[mn + "_sd"].append(0)
                    output_dict[mn + "_mad_to_sd_ratio"].append(0)
                    for i in range(1, 4):
                        output_dict[f"{mn}_{i}MAD"].append(0)
                else:
                    median = cluster_cells[mn + "_median"][0]
                    mad = cluster_cells[mn + "_mad"][0]
                    sd = np.std(cluster_cells[mn])

                    output_dict[mn + "_median"].append(median)
                    output_dict[mn + "_mad"].append(mad)
                    output_dict[mn + "_sd"].append(sd)
                    output_dict[mn + "_mad_to_sd_ratio"].append(mad / sd)
                    for i in range(1, 4):
                        if mn == "percent_mito":
                            co = median + i * mad
                            # if co > max(cluster_cells[mn]):
                            #     co = None
                        else:
                            co = median - i * mad
                            # if co < min(cluster_cells[mn]):
                            #     co = None
                        output_dict[f"{mn}_{i}MAD"].append(co)

    output_path.mkdir(parents=True, exist_ok=True)
    with open(output_path / "filtered_cells_count.csv", "w") as file:
        df = pd.DataFrame.from_dict(output_dict)
        file.write(df.to_csv())
    with open(output_path / "initial_clustering.csv", "w") as file:
        file.write(cells.to_csv())

    # launch R plot script
    print(r"& 'C:\Program Files\R\R-4.1.3\bin\Rscript.exe' {} {}-{}".format(FILE_PATH / "adaptive_mad_threshold_plot.R", project, tissue))
