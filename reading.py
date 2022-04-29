import os

import pandas as pd
import pegasusio as io

from config.config import DATA_DIR


# function that parses projects.csv and returns relevant info
def get_project_info(project=None, task_id=None, tissue=None):
    try:
        projects = pd.read_csv("config/read_info/projects.csv")
    except FileNotFoundError:
        projects = pd.read_csv("../config/read_info/projects.csv")
    if project is None:
        return projects
    assert project in set(projects['project'])  # check if project exists in project list
    project_tissues = projects[projects['project'] == project]  # subset tissues for the project
    if task_id:  # get info based on task_id
        assert task_id < len(project_tissues['project'])
        tissue = list(project_tissues['tissue'])[task_id]
        is_human = list(project_tissues['is_human'])[task_id]
        annotations = list(project_tissues['annotations'])[task_id]
        return tissue, is_human, annotations
    elif tissue:  # get info based on tissue
        assert tissue in set(project_tissues['tissue'])
        tissue_info = project_tissues[project_tissues['tissue'] == tissue]
        is_human = list(tissue_info['is_human'])[0]
        annotations = list(tissue_info['annotations'])[0]
        return is_human, annotations
    else:  # if no tissue or task_id provided return pandas project with project info
        return project_tissues


def read_tissue(project, tissue, annotations="Unknown"):  # function that reads and aggregates one tissue
    try:
        dataset_list_path = "config/read_info/{}/{}.csv".format(project, tissue)  # path to tissue read info
        assert os.path.isfile(dataset_list_path)  # check if tissue read info exists
    except AssertionError:
        dataset_list_path = "../config/read_info/{}/{}.csv".format(project, tissue)  # path to tissue read info
        assert os.path.isfile(dataset_list_path)  # check if tissue read info exists
    read_info_filename = "read_info_{}_{}.csv".format(project, tissue)  # filename of a current read info copy
    dataset_list = pd.read_csv(dataset_list_path)
    if annotations == "Absolute":  # if data is stored outside DATA_DIR and csv had a direct path
        pass
    else:
        dataset_list['Location'] = [DATA_DIR + t for t in dataset_list['Location']]  # update location with relevant directory prefix
    with open(read_info_filename, "w") as fout:  # write modified csv
        fout.write(dataset_list.to_csv())
    adata = io.aggregate_matrices(read_info_filename)
    os.remove(read_info_filename)  # remove current read info copy

    # add annotations to adata
    if annotations != 'Unknown' and annotations != "Absolute" and os.path.isfile(DATA_DIR + annotations):
        ann_df = pd.read_csv(DATA_DIR + annotations)
        annotations_cell_type = ann_df["annotations"]
        annotations_cell_type.index = ann_df["barcodekey"]
        annotations_cell_type = annotations_cell_type.reindex(adata.obs.index)
        annotations_cell_type = annotations_cell_type.fillna("Unknown")
        adata.obs["annotations"] = annotations_cell_type
    else:
        adata.obs['annotations'] = 'Unknown'
    return adata
