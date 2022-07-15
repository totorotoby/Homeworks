#!/bin/bash

for (( n = 2 ; n <= 24 ; n++ ))
do
    export OMP_NUM_THREADS=$n
    ./prefix 1048576 1
done
