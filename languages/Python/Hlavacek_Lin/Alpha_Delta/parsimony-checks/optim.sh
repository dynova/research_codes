#!/bin/bash
#SBATCH --chdir=/home/dynova/n3
#SBATCH --time=02:00:00
#SBATCH --job-name=nm
#SBATCH --array=0-1499
#SBATCH --output=nm-%a.out
#SBATCH --error=nm-%a.err
#SBATCH --mail-type=END
#SBATCH --mail-user=abhishek.mallela@gmail.com
python3 optim.py $SLURM_ARRAY_TASK_ID