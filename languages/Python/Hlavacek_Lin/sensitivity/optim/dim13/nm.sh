#!/bin/bash
#SBATCH --chdir=/home/dynova/optim/dim13
#SBATCH --time=15:00:00
#SBATCH --job-name=nm
#SBATCH --array=0-99
#SBATCH --output=nm-%a.out
#SBATCH --error=nm-%a.err
#SBATCH --mail-type=END
#SBATCH --mail-user=abhishek.mallela@gmail.com
python3 nm.py $SLURM_ARRAY_TASK_ID