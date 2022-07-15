#!/bin/bash

for (( n = 2 ; n <= 1024 ; n = n * 2))
do
    ./prefix $n 1
done
