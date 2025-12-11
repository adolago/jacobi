#!/bin/bash

[ -d results ] || mkdir results

tol=1e-6
max_iters=1000000

for ((i=3;i<5;i++));do
   nx=$((2**i+1))
   odir=`printf "results/NX%04d" ${nx}`
   [ -d ${odir} ] || mkdir ${odir}
   ./JACOBI2D ${nx} ${tol} ${max_iters} ${odir}
done
