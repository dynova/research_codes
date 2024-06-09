#!/bin/bash
#PBS -N PDT2
#PBS -l nodes=1:ppn=1,mem=5g,walltime=48:00:00
#PBS -m abe

g++ PDT2.cpp mtrand.cpp -o PDT2
./PDT2 1000000 100 var_PDS2PDT2_unsat.txt PDT2_Unsat.txt