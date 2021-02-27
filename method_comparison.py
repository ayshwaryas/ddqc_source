import subprocess
import sys

from config.config import *
from filtering import filter_cells
from reading import read_tissue, get_project_info
from utils import cluster_data, safe_mkdir, add_cd_scores, write_markers, save_to_csv


# function that creates all the relevant directories
def create_dirs(project, tissue, res, method, param):
    task_directory = "{}-{}-{}".format(res, param, method)  # name of the directory for this task
    task_name = tissue + "-" + task_directory  # task name to put on plots
    results_dir = OUTPUT_DIR + project + "/" + tissue + "/" + task_directory + "/"  # directory for saving output

    # create all the directories, if they dont exist
    print("Creating Output Directories")
    safe_mkdir(OUTPUT_DIR)
    safe_mkdir(OUTPUT_DIR + project + "/")
    safe_mkdir(OUTPUT_DIR + project + "/" + tissue)
    safe_mkdir(results_dir)

    return task_directory, task_name, results_dir


def mc_main(project, task_id, tissue=None):
    if tissue is None:
        tissue, is_human, annotations = get_project_info(project, task_id=task_id // MC_TASKS_PER_TISSUE)
    else:
        is_human, annotations = get_project_info(project, tissue=tissue)

    adata = read_tissue(project, tissue)

    method, param = MC_METHODS[task_id % MC_TASKS_PER_TISSUE]
    print("Method comparison task_id:{} - tissue:{}, res:{}, method:{}, param:{}, project:{}".format(task_id, tissue,
                                                                                                     resolution, method,
                                                                                                     param, project,
                                                                                                     do_counts,
                                                                                                     do_genes, do_mito,
                                                                                                     do_ribo))

    task_directory, task_name, results_dir = create_dirs(project, tissue, resolution, method, param)

    # filtering
    mito_prefix = MITO_PREFIXES["human"] if is_human else MITO_PREFIXES["mouse"]
    ribo_prefix = RIBO_PREFIXES["human"] if is_human else RIBO_PREFIXES["mouse"]
    adata = filter_cells(adata, resolution, method, param, basic_n_genes=basic_genes_filter,
                         basic_percent_mito=basic_mito_filter, mito_prefix=mito_prefix, ribo_prefix=ribo_prefix,
                         do_counts=do_counts, do_genes=do_genes, do_mito=do_mito, do_ribo=do_ribo)

    adata, marker_dict = cluster_data(adata, compute_markers=True, compute_reductions=True, resolution=resolution)
    adata = add_cd_scores(adata, is_human)  # add cell death scores

    # write the results
    print("Writing results")
    write_markers(marker_dict, results_dir)
    save_to_csv(adata, results_dir)

    # launch R plot script
    print(subprocess.check_output("Rscript r_plots.R {} {} {} {} {}".format(project, tissue, resolution, method, param),
                                  shell=True).decode('UTF-8'))


if __name__ == '__main__':
    if local:  # for debug outside of cluster
        proj = input("Project: ").strip()
        t_id = int(input("Task ID (-1 to input specific tissue): ").strip())
        if t_id == -1:
            tiss = input("Tissue: ").strip()
            mc_main(proj, 0, tissue=tiss)
        else:
            mc_main(proj, t_id)
    else:  # project and task id are provided as commandline args
        proj = sys.argv[1]
        t_id = int(sys.argv[2]) - 1
        mc_main(proj, t_id)