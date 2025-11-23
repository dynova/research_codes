#!/bin/bash
#PBS -N PDS2
#PBS -l nodes=1:ppn=1,mem=5g,walltime=48:00:00
#PBS -m abe

g++ PDS2.cpp mtrand.cpp -o PDS2
./PDS2 1000000 100 var_PDS2PDT2_unsat.txt PDS2_Unsat.txt