#!/bin/bash -l
#SBATCH --chdir=/vast/home/amallela/codes/EID/CFD/LSODA/withoutJacobian
#SBATCH --nodes=1
#SBATCH --constraint="cpu_family:haswell" 
#SBATCH --qos=normal
#SBATCH --mem=5G
#SBATCH --ntasks-per-node=8
#SBATCH --time=10:00:00
#SBATCH --job-name=cfd
#SBATCH --exclusive
#SBATCH --mail-type=END
#SBATCH --export=OPENBLAS_NUM_THREADS=1
#SBATCH --mail-user=amallela@lanl.gov
module load miniconda3
source activate jupyter-env
python3 cfd.py
lscpu > cpuInfo_cfd.txt
cat /proc/meminfo > memoryInfo_cfd.txt