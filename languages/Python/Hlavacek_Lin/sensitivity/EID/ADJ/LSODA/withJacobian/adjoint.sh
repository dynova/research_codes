#!/bin/bash -l
#SBATCH --chdir=/vast/home/amallela/codes/EID/ADJ/LSODA/withJacobian
#SBATCH --nodes=1
#SBATCH --constraint="cpu_family:haswell" 
#SBATCH --qos=normal
#SBATCH --mem=5G
#SBATCH --ntasks-per-node=8
#SBATCH --time=10:00:00
#SBATCH --job-name=adjoint
#SBATCH --exclusive
#SBATCH --mail-type=END
#SBATCH --export=OPENBLAS_NUM_THREADS=1
#SBATCH --mail-user=amallela@lanl.gov
module load miniconda3
source activate jupyter-env
python3 adjoint.py
lscpu > cpuInfo_adjoint.txt
cat /proc/meminfo > memoryInfo_adjoint.txt