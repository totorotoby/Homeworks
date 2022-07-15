#!/bin/bash

for (( n = 1 ; n <= 12 ; n++ ))
do
    export OMP_NUM_THREADS=$n
    ./spmv cant/cant.mtx cant/b.mtx cant/myanscant.matrix
done
