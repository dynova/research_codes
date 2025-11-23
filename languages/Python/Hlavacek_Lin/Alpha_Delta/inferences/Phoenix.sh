#!/bin/bash -l
#SBATCH --chdir=/vast/home/amallela/test
#SBATCH --time=01:00:00
#SBATCH --qos=debug
#SBATCH --job-name=P
#SBATCH --array=0-31
#SBATCH --mail-type=END
#SBATCH --mail-user=amallela@lanl.gov
module load miniconda3
source activate jupyter-env
python3 Phoenix.py $SLURM_ARRAY_TASK_ID