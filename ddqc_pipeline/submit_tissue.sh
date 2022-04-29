#!/bin/bash

#$ -l h_vmem=15G
#$ -l h_rt=48:00:00
#$ -P regevlab
# Cores
#$ -pe smp 6
#$ -R y
#$ -binding linear:6

#$ -o /ahg/regevdata/projects/scqc/code/logs/

source /broad/software/scripts/useuse
use .python-3.8.3
use R-3.5

source /home/unix/malperov/pegasusenv/bin/activate

cd /home/unix/malperov/Primes2019/ddqc_pipeline/
python submission/run_tissue.py $1 $2 $3 $4