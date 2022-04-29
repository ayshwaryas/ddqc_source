import os
import sys

from reading import get_project_info


method = sys.argv[1]
run_analysis = sys.argv[2]
projects = get_project_info()

for i in projects.iterrows():
    project = i[1]["project"]
    tissue = i[1]["tissue"]
    os.system(f"qsub submit_tissue.sh {project} {tissue} {method} {run_analysis}")
