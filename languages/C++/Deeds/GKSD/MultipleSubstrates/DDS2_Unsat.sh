#!/bin/bash
#PBS -N DDS2
#PBS -l nodes=1:ppn=1,mem=1g,walltime=16:00:00
#PBS -m abe

g++ DDS2.cpp mtrand.cpp -o DDS2
./DDS2 1000000 100 var_DDS2DDT2_unsat.txt DDS22.txt