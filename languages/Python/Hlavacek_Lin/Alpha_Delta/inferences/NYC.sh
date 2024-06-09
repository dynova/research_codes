#!/bin/bash
#SBATCH --chdir=/home/dynova/test
#SBATCH --time=04:00:00
#SBATCH --job-name=N
#SBATCH --array=0-24
#SBATCH --mail-type=END
#SBATCH --mail-user=abhishek.mallela@gmail.com
python3 NYC.py $SLURM_ARRAY_TASK_ID