1#!/bin/bash

for (( n = 1 ; n <= 12 ; n++ ))
do
    export OMP_NUM_THREADS=$n
    ./cg cant/cant.mtx cant/b.mtx cant/myres.matrix 13000 .000001 1
done
