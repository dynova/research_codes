#!/bin/bash
#SBATCH --chdir=/home/dynova/optim/dim7
#SBATCH --time=15:00:00
#SBATCH --job-name=gd
#SBATCH --array=0-99
#SBATCH --output=gd-%a.out
#SBATCH --error=gd-%a.err
#SBATCH --mail-type=END
#SBATCH --mail-user=abhishek.mallela@gmail.com
python3 gd.py $SLURM_ARRAY_TASK_ID