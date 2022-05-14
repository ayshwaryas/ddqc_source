import pandas as pd


path = "Z:\\revisions\\revisions2\\cell_classification"
tissue = "krasnow_lung"

cells = pd.read_csv(f"{path}\\{tissue}\\classification_mad_results.csv", index_col=0)
ddqc_cells = cells[cells.ddqc_cluster_passed_qc]
cell_typist_cells = cells[cells.cell_typist_passed_qc]


# for cluster in sorted(pd.unique(ddqc_cells.ddqc_cluster)):
#     cluster_cells = ddqc_cells[ddqc_cells.ddqc_cluster == cluster]
#     n_cells = len(cluster_cells.index)
#
#     for i in range(3):
#         ct_cell_type = cluster_cells.cell_typist.value_counts().index[i]
#         ct_n_cells = cluster_cells.cell_typist.value_counts()[ct_cell_type]
#         print(f"{i + 1}. ddqc cluster {cluster}: {ct_n_cells} out of {n_cells} ({round(ct_n_cells / n_cells * 100, 3)}%) are {ct_cell_type}")
#
#     n_cells_passed_ct = len(cluster_cells[cluster_cells.cell_typist_passed_qc].index)
#     print(f"ddqc cluster {cluster}: {n_cells_passed_ct} out of {n_cells} ({round(n_cells_passed_ct / n_cells * 100, 3)}%) passed Cell Typist")
#     print()


for cell_type in sorted(pd.unique(cell_typist_cells.cell_typist)):
    cell_type_cells = cell_typist_cells[cell_typist_cells.cell_typist == cell_type]
    n_cells = len(cell_type_cells.index)

    for i in range(min(3, len(cell_type_cells.ddqc_cluster.value_counts().index))):
        ddqc_cluster = cell_type_cells.ddqc_cluster.value_counts().index[i]
        ddqc_n_cells = cell_type_cells.ddqc_cluster.value_counts()[ddqc_cluster]
        print(f"{i + 1}. Cell Typist {cell_type}: {ddqc_n_cells} out of {n_cells} ({round(ddqc_n_cells / n_cells * 100, 3)}%) are ddqc cluster {ddqc_cluster}")

    n_cells_passed_ddqc = len(cell_type_cells[cell_type_cells.ddqc_cluster_passed_qc].index)
    print(f"Cell Typist {cell_type}: {n_cells_passed_ddqc} out of {n_cells} ({round(n_cells_passed_ddqc / n_cells * 100, 3)}%) passed ddqc")
    print()
