from mc import main as mc_main
from mc_plot import main as mc_plot_main
import traceback

project = "mc_heart_circulation"
for task_id in range(1 * 3):
    try:
        mc_main(project, task_id)
    except Exception as e:
        print("========= FAIL TASK_ID %s =========" % task_id)
        print(traceback.format_exc())
for task_id in range(1 * 1):
    try:
        mc_plot_main(project, task_id)
    except Exception as e:
        print("========= FAIL TASK_ID %s =========" % task_id)
        print(traceback.format_exc())