import sys

from method_comparison import mc_main
from joint_clustering_old import joint_main
from config.config import MC_TASKS_PER_TISSUE

project = "tabula_muris"
tissues = ["Heart_and_Aorta", "Lung"]
method = 2
for tissue in tissues:
    mc_main(project, task_id=method, tissue=tissue, method=method)
