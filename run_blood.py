from method_comparison import mc_main
from joint_clustering import joint_main
from config.config import MC_TASKS_PER_TISSUE

project = "human_other"
tissue = "Blood"
for method in range(MC_TASKS_PER_TISSUE):
    mc_main(project, task_id=method, tissue=tissue)
joint_main(project, 0, tissue=tissue)