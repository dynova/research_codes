#!/bin/bash
#SBATCH --chdir=/home/dynova/optim/dim7
#SBATCH --time=15:00:00
#SBATCH --job-name=tc
#SBATCH --array=0-99
#SBATCH --output=tc-%a.out
#SBATCH --error=tc-%a.err
#SBATCH --mail-type=END
#SBATCH --mail-user=abhishek.mallela@gmail.com
python3 tc.py $SLURM_ARRAY_TASK_ID