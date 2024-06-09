#!/bin/bash -l
#SBATCH --chdir=/vast/home/amallela/codes/EID/AD/EULER
#SBATCH --nodes=1
#SBATCH --constraint="cpu_family:haswell" 
#SBATCH --qos=debug
#SBATCH --mem=5G
#SBATCH --ntasks-per-node=8
#SBATCH --time=04:00:00
#SBATCH --job-name=euler
#SBATCH --exclusive
#SBATCH --mail-type=END
#SBATCH --export=OPENBLAS_NUM_THREADS=1
#SBATCH --mail-user=amallela@lanl.gov
module load miniconda3
source activate jupyter-env
python3 euler.py
lscpu > cpuInfo_euler.txt
cat /proc/meminfo > memoryInfo_euler.txt