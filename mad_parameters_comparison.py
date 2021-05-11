import sys

from method_comparison import mc_main
from joint_clustering_old import joint_main
from config.config import MC_TASKS_PER_TISSUE

project = "tabula_muris"
method = 2
for tissue in ("Heart_and_Aorta", "Lung"):
    for param in (1.5, 2.5):
        mc_main(project, task_id=method, tissue=tissue, param=param)
