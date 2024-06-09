#!/bin/bash -l
#SBATCH --chdir=/vast/home/amallela/test
#SBATCH --time=04:00:00
#SBATCH --qos=debug
#SBATCH --job-name=H
#SBATCH --array=0-15
#SBATCH --mail-type=END
#SBATCH --mail-user=amallela@lanl.gov
module load miniconda3
source activate jupyter-env
python3 Houston.py $SLURM_ARRAY_TASK_ID