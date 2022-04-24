import subprocess

from config.config import *
from filtering import filter_cells
from reading import read_tissue, get_project_info
from utils import cluster_data, safe_mkdir, add_cd_scores, marker_dict_to_df, save_to_csv, assign_cell_types


# function that creates all the relevant directories
def create_dirs(project, tissue, clustering_params, res, method, param):
    task_directory = "{}-{}-{}-{}-{}".format(res, clustering_params[1], clustering_params[2], method,
                                             param)  # name of the directory for this task
    task_name = tissue + "-" + task_directory  # task name to put on plots
    results_dir = OUTPUT_DIR + project + "/" + tissue + "/" + task_directory + "/"  # directory for saving output

    # create all the directories, if they dont exist
    print("Creating Output Directories")
    safe_mkdir(OUTPUT_DIR)
    safe_mkdir(OUTPUT_DIR + project + "/")
    safe_mkdir(OUTPUT_DIR + project + "/" + tissue)
    safe_mkdir(results_dir)

    return task_directory, task_name, results_dir


def mc_main(project, tissue, clustering_params, task_id):
    is_human, annotations = get_project_info(project, tissue=tissue)

    adata = read_tissue(project, tissue, annotations)
    method, param = MC_METHODS[task_id % MC_TASKS_PER_TISSUE]
    resolution = clustering_params[0]
    print("Clustering comparison - clustering_params:{}, tissue:{}, res:{}, method:{}, param:{}, project:{}".format(
        clustering_params, tissue,
        resolution, method,
        param, project))

    task_directory, task_name, results_dir = create_dirs(project, tissue, clustering_params, resolution, method, param)

    # filtering
    mito_prefix = MITO_PREFIXES["human"] if is_human else MITO_PREFIXES["mouse"]
    ribo_prefix = RIBO_PREFIXES["human"] if is_human else RIBO_PREFIXES["mouse"]
    adata = filter_cells(adata, resolution, method, param, basic_n_genes=basic_genes_filter,
                         basic_percent_mito=basic_mito_filter, mito_prefix=mito_prefix, ribo_prefix=ribo_prefix,
                         do_counts=do_counts, do_genes=do_genes, do_mito=do_mito, do_ribo=do_ribo,
                         record_path=results_dir, n_components=clustering_params[1], k=clustering_params[2])

    adata, marker_dict = cluster_data(adata, compute_markers=True, compute_reductions=True, resolution=resolution,
                                      n_components=clustering_params[1], k=clustering_params[2])
    adata = add_cd_scores(adata, is_human)  # add cell death scores
    markers = marker_dict_to_df(marker_dict)
    clusters = assign_cell_types(adata, markers, tissue)

    # write the results
    print("Writing results")
    save_to_csv(adata, results_dir)
    with open(results_dir + "!clusters.csv", "w") as fout:
        fout.write(clusters.to_csv())
    with open(results_dir + "!markers.csv", "w") as fout:
        fout.write(markers.to_csv())

    # launch R plot script
    print(
        subprocess.check_output("Rscript plotting.R {} {} {}".format(task_name, results_dir, "mc"),
                                shell=True).decode('UTF-8'))


if __name__ == '__main__':
    if local:  # for debug outside of cluster
        proj = "tabula_muris"
        tiss = "Lung"
        t_id = 2
        for clust_param in [[1.4, 50, 20]]: #(
        #         [1, 50, 20], [0.5, 50, 20], [1.4, 50, 10], [1, 50, 10], [0.5, 50, 10], [1.4, 20, 20], [1, 20, 20],
        #         [0.5, 20, 20], [1.4, 20, 10], [1, 20, 10], [0.5, 20, 10],):
            mc_main(proj, tiss, clust_param, int(t_id))

    # else:  # project and task id are provided as commandline args
    #     proj = sys.argv[1]
    #     t_id = int(sys.argv[2]) - 1
    #     mc_main(proj, t_id)
