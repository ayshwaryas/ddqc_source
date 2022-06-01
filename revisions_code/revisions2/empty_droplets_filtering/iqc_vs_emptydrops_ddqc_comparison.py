import pandas as pd


path = "Z:\\output_pg\\tabula_muris\\Heart_and_Aorta\\1.4-mad-2"

iqc_cells = pd.read_csv(f"{path}\\!cells_initial.csv", index_col=0)
emptydrops_cells = pd.read_csv(f"{path}\\!cells_initial_emptydrops.csv", index_col=0)
common_barcodes = list(set(iqc_cells.index).intersection(set(emptydrops_cells.index)))

iqc_cells = iqc_cells[[t in common_barcodes for t in iqc_cells.index]]
iqc_cells = iqc_cells[iqc_cells.passed_qc]
emptydrops_cells = emptydrops_cells[[t in common_barcodes for t in emptydrops_cells.index]]
emptydrops_cells = emptydrops_cells[emptydrops_cells.passed_qc]

print("Total cells IQC after ddqc:", len(iqc_cells.index))
print("Total cells emptydrops after ddqc:", len(emptydrops_cells.index))
print("Cells in common:", len(set(iqc_cells.index).intersection(set(emptydrops_cells.index))))
print("Cells unique to IQC:", len(set(iqc_cells.index).difference(set(emptydrops_cells.index))))
print("Cells unique to emptydrops:", len(set(emptydrops_cells.index).difference(set(iqc_cells.index))))
print("\n\n")

for cluster in sorted(pd.unique(iqc_cells.cluster_labels)):
    cluster_cells = iqc_cells[iqc_cells.cluster_labels == cluster]
    n_cells = len(cluster_cells.index)
    emptydrops_cluster_cells = emptydrops_cells[[t in cluster_cells.index for t in emptydrops_cells.index]]

    for i in range(min(3, len(emptydrops_cluster_cells.cluster_labels.value_counts().index))):
        emptydrops_cl = emptydrops_cluster_cells.cluster_labels.value_counts().index[i]
        emptydrops_n_cells = emptydrops_cluster_cells.cluster_labels.value_counts()[emptydrops_cl]
        print(f"{i + 1}. iqc cluster {cluster}: {emptydrops_n_cells} out of {n_cells} "
              f"({round(emptydrops_n_cells / n_cells * 100, 3)}%) are from emptydrops cluster {emptydrops_cl}")

    n_cells_passed_emptydrops = len(emptydrops_cluster_cells.index)
    print(f"iqc cluster {cluster}: {n_cells_passed_emptydrops} out of {n_cells} "
          f"({round(n_cells_passed_emptydrops / n_cells * 100, 3)}%) passed emptydrops ddqc")
    print()

print("\n====================================\n")

for cluster in sorted(pd.unique(emptydrops_cells.cluster_labels)):
    cluster_cells = emptydrops_cells[emptydrops_cells.cluster_labels == cluster]
    n_cells = len(cluster_cells.index)
    iqc_cluster_cells = iqc_cells[[t in cluster_cells.index for t in iqc_cells.index]]

    for i in range(min(3, len(iqc_cluster_cells.cluster_labels.value_counts().index))):
        iqc_cl = iqc_cluster_cells.cluster_labels.value_counts().index[i]
        iqc_n_cells = iqc_cluster_cells.cluster_labels.value_counts()[iqc_cl]
        print(f"{i + 1}. emptydrops cluster {cluster}: {iqc_n_cells} out of {n_cells} "
              f"({round(iqc_n_cells / n_cells * 100, 3)}%) are from iqc cluster {iqc_cl}")

    n_cells_passed_iqc = len(iqc_cluster_cells.index)
    print(f"emptydrops cluster {cluster}: {n_cells_passed_iqc} out of {n_cells} "
          f"({round(n_cells_passed_iqc / n_cells * 100, 3)}%) passed iqc ddqc")
    print()
