#!/bin/bash
#PBS -N PP2
#PBS -q long
#PBS -l nodes=1:ppn=1,mem=20g,walltime=672:00:00
#PBS -m abe

g++ PP2.cpp mtrand.cpp -o PP2
./PP2 1000000 100 var_PP2_unsat.txt PP2_Unsat.txt