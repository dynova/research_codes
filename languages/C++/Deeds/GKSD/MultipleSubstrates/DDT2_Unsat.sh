#!/bin/bash
#PBS -N DDT2
#PBS -l nodes=1:ppn=1,mem=1g,walltime=16:00:00
#PBS -m abe

g++ DDT2.cpp mtrand.cpp -o DDT2
./DDT2 1000000 100 var_DDS2DDT2_unsat.txt DDT22.txt