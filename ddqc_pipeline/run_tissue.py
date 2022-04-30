import sys

from method_comparison import mc_main
from joint_clustering_old import joint_main
from config.config import MC_TASKS_PER_TISSUE

project = sys.argv[1]
tissue = sys.argv[2]
method_id = int(sys.argv[3])
run_analysis = sys.argv[4] == "True"

if method_id == -1:
    for method in range(MC_TASKS_PER_TISSUE):
        mc_main(project, task_id=method, tissue=tissue, run_analysis=run_analysis)
    joint_main(project, 0, tissue=tissue)
else:
    mc_main(project, task_id=method_id, tissue=tissue, run_analysis=run_analysis)
