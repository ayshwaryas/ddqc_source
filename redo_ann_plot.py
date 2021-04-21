import subprocess
from collections import Counter

import numpy as np
import pandas as pd

from config.config import *
from joint_clustering_old import create_joint_dir
from method_comparison import create_dirs
from reading import get_project_info
from utils import cell_type_replacements, format_markers


def assign_cell_types(adata, markers, tissue, min_sg=3):
    print("Assigning annotations and cell types to clusters")
    genes = pd.read_csv(DATA_DIR + "markers.tsv", sep="\t")  # read panglaoDB genes to cell types table
    colnames = ["cluster", "annotation", "annotation2", "%annotation", "%annotation2", "cell_type", "cell_type2",
                "n_cells", "genes_mean", "genes_median", "mito_mean", "mito_median", "ribo_mean", "ribo_median",
                "markers", "score_genes"]  # colnames for cluster info csv
    clusters_dict = dict()  # create a dict with colnames, that will be later converted to a df
    for c in colnames:
        clusters_dict[c] = []

    for cl in range(1, max(list(set(adata.louvain_labels))) + 1):  # iterate though all clusters
        cluster = adata[adata.louvain_labels == cl]  # subset adata based on cluster number
        cluster_markers = markers[markers.cluster == cl]  # subset markers for current cluster

        possible_cell_types = dict()  # dict with possible cell types and scores (sum of Log2FC for each marker)
        score_genes = dict()  # dict with possible cell types and score_genes for them
        for index, marker_info in cluster_markers.iterrows():  # iterate through markers
            marker = index.upper()
            gene_info = genes[genes["official gene symbol"] == marker]  # find marker in PanglaoDB table
            for gene, info in gene_info.iterrows():  # iterate through  cell types for the marker
                ct = cell_type_replacements(info["cell type"], tissue)
                if ct not in possible_cell_types:  # if we didn't encounter this cell type yet, add it to dicts
                    possible_cell_types[ct] = 0
                    score_genes[ct] = []
                possible_cell_types[ct] += marker_info["log2FC"]  # increase score by Log2FC of current marker
                score_genes[ct].append(marker)  # add marker to the list of score genes for current cell type

        called = False  # whether the cell type for the cluster was called
        if len(possible_cell_types) > 1:  # more than 2 cell types
            # sort cell types by score
            srt_cell_types = sorted(list(possible_cell_types.keys()), key=lambda x: -possible_cell_types[x])
            sg = score_genes[srt_cell_types[0]]
            if len(sg) >= min_sg:  # if number of score genes satisfies the minimum requirement call the cell type
                clusters_dict["cell_type"].append(srt_cell_types[0])
                clusters_dict["cell_type2"].append(srt_cell_types[1])
                clusters_dict["score_genes"].append(format_markers(sg))
                called = True
        elif len(possible_cell_types) == 1:  # only one cell type
            ct = list(possible_cell_types.keys())[0]
            sg = score_genes[ct]
            if len(sg) >= min_sg:  # if number of score genes satisfies the minimum requirement call the cell type
                clusters_dict["cell_type"].append(ct)
                clusters_dict["cell_type2"].append("")  # second common cell type doesn't exist
                clusters_dict["score_genes"].append(format_markers(sg))
                called = True
        if not called:  # if cell type wasn't called, mark it Unknown
            clusters_dict["cell_type"].append("Unknown")
            clusters_dict["cell_type2"].append("")
            clusters_dict["score_genes"].append([])

        # calculate first and second most frequent annotation and percentage of cells that have them
        cl_anno = list(cluster.annotations)
        most_common = Counter(cl_anno).most_common()
        clusters_dict["annotation"].append(most_common[0][0])
        clusters_dict["%annotation"].append(round(most_common[0][1] / len(cl_anno), 3))
        if len(most_common) > 1:  # more than one annotation
            clusters_dict["annotation2"].append(most_common[1][0])
            clusters_dict["%annotation2"].append(round(most_common[1][1] / len(cl_anno), 3))
        else:  # only one annotation
            clusters_dict["annotation2"].append("")
            clusters_dict["%annotation2"].append(0)

        # record cluster statistics and markers
        clusters_dict["cluster"].append(cl)
        clusters_dict["n_cells"].append(len(cl_anno))
        clusters_dict["genes_mean"].append(round(float(np.mean(cluster.n_genes)), 3))
        clusters_dict["genes_median"].append(round(float(np.median(cluster.n_genes)), 3))
        clusters_dict["mito_mean"].append(round(float(np.mean(cluster.percent_mito)), 3))
        clusters_dict["mito_median"].append(round(float(np.median(cluster.percent_mito)), 3))
        clusters_dict["ribo_mean"].append(round(float(np.mean(cluster.percent_ribo)), 3))
        clusters_dict["ribo_median"].append(round(float(np.median(cluster.percent_ribo)), 3))

        clusters_dict["markers"].append(format_markers(list(cluster_markers.index)))

        print("Cluster {} finished".format(cl))

    return pd.DataFrame.from_dict(clusters_dict)  # convert dict to pandas df


def remake(project, task_id, tissue=None):
    if tissue is None:
        tissue = get_project_info(project, task_id=task_id // MC_TASKS_PER_TISSUE)[0]

    if task_id < 3:
        method, param = MC_METHODS[task_id % MC_TASKS_PER_TISSUE]
        print(
            "Method comparison task_id:{} - tissue:{}, res:{}, method:{}, param:{}, project:{}".format(task_id, tissue,
                                                                                                       resolution,
                                                                                                       method, param,
                                                                                                       project))
        task_directory, task_name, results_dir = create_dirs(project, tissue, resolution, method, param)
    else:
        print(
            "Joint clustering task_id:{} - tissue:{}, res:{}, project:{}".format(task_id, tissue, resolution, project))
        task_directory, task_name, results_dir = create_joint_dir(project, tissue, resolution)

    source_dir = OUTPUT_DIR.split("output_pg")[0] + "output_copies/output_pg04-20-21/" + results_dir.split(OUTPUT_DIR)[
        1]
    adata = pd.read_csv(source_dir + "!cells.csv")
    markers = pd.read_csv(source_dir + "!markers.csv")
    markers.index = markers.feature
    if not any(markers.t_qval < 0.05):
        markers = markers[markers.t_qval < 0.1]
        print("Increased q_val threshold to 0.1")
    else:
        markers = markers[markers.t_qval < 0.05]

    clusters = assign_cell_types(adata, markers, tissue)

    with open(results_dir + "!cells.csv", "w") as fout:
        fout.write(adata.to_csv())
    with open(results_dir + "!clusters.csv", "w") as fout:
        fout.write(clusters.to_csv())
    with open(results_dir + "!markers.csv", "w") as fout:
        fout.write(markers.to_csv())

    # launch R plot script
    if task_id < 3:
        print(
            subprocess.check_output("Rscript plotting.R {} {} {}".format(task_name, results_dir, "mc"),
                                    shell=True).decode('UTF-8'))
    else:
        print(
            subprocess.check_output("Rscript plotting.R {} {} {}".format(task_name, results_dir, "joint"),
                                    shell=True).decode('UTF-8'))


for task in ("human_other/Adipose",
             "human_other/Kidney2",
             "human_other/Liver",
             "human_other/Krasnow_Lung",
             "PanglaoDB/Bone_Marrow",
             "PanglaoDB/Mammary_Gland",
             "PanglaoDB/Pancreatic_Islets",
             "PanglaoDB/Substantia_Nigra",
             "PanglaoDB/Testis",
             "tabula_muris_smartseq2/Bone_Marrow",
             "tabula_muris_smartseq2/Cerebellum",
             "tabula_muris_smartseq2/Colon",
             "tabula_muris_smartseq2/Heart_and_Aorta",
             "tabula_muris/Bladder",
             "tabula_muris/Heart_and_Aorta",
             "tabula_muris/Lung",
             "tabula_muris/Mammary_Gland",
             "tabula_muris/Tongue",
             "tabula_muris/Trachea",
             "human_other/Heart_Circulation"):

    proj, tiss = task.split("/")
    for method in range(4):
        remake(proj, method, tiss)
