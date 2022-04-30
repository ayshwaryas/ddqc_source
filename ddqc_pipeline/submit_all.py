import os
import pathlib
import sys

from config.config import DATA_DIR
from reading import get_project_info

PATH = pathlib.Path(DATA_DIR).parent / "submit_info"
(PATH / "logs").mkdir()
(PATH / "scripts").mkdir()

COMMAND = """#!/bin/bash

#$ -l h_vmem=15G
#$ -l h_rt=48:00:00
# Cores
#$ -pe smp 6
#$ -R y
#$ -binding linear:6

#$ -o {path}/logs/
#$ -N {project}-{tissue}-{method}

source /broad/software/scripts/useuse
use .python-3.8.3
use R-3.5

source /home/unix/malperov/pegasusenv/bin/activate
cd /home/unix/malperov/Primes-2019/ddqc_pipeline/
python submission/run_tissue.py {project} {tissue} {method} {run_analysis}"""


method = sys.argv[1]
run_analysis = sys.argv[2]
projects = get_project_info()

for i in projects.iterrows():
    project = i[1]["project"]
    tissue = i[1]["tissue"]
    cmd = COMMAND.format(path=PATH, project=project, tissue=tissue, method=method, run_analysis=run_analysis)
    if project == "tabula_muris" and tissue == "Lung":
        file = PATH / f"scripts/{project}_{tissue}_{method}.sh"
        file.write_text(cmd)
        os.system(f"qsub {file.absolute()}")
