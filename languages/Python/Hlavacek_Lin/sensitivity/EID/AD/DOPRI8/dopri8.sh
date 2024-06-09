#!/bin/bash -l
#SBATCH --chdir=/vast/home/amallela/codes/EID/AD/DOPRI8
#SBATCH --nodes=1
#SBATCH --constraint="cpu_family:haswell" 
#SBATCH --qos=debug
#SBATCH --mem=5G
#SBATCH --ntasks-per-node=8
#SBATCH --time=04:00:00
#SBATCH --job-name=dopri8
#SBATCH --exclusive
#SBATCH --mail-type=END
#SBATCH --export=OPENBLAS_NUM_THREADS=1
#SBATCH --mail-user=amallela@lanl.gov
module load miniconda3
source activate jupyter-env
python3 dopri8.py
lscpu > cpuInfo_dopri8.txt
cat /proc/meminfo > memoryInfo_dopri8.txt