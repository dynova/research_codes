#!/bin/bash
#PBS -N DP2
#PBS -q long
#PBS -l nodes=1:ppn=1,mem=20g,walltime=672:00:00
#PBS -m abe

g++ DP2.cpp mtrand.cpp -o DP2
./DP2 1000000 100 var_DP2_unsat.txt DP2_Unsat.txt