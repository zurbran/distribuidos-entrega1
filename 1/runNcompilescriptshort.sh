#!/bin/bash

cd 3
cd openmp
make
make sequential

for size in "200000000" "400000000" "800000000"
do
    for threads in 2 4
    do
        for runs in 1 2 3 4
        do
            echo -e "./2_openmp $size $threads"
            ./3_openmp $size $threads
        done
    done
done

for size in "200000000" "400000000" "800000000"
do
    for runs in 1 2 3 4
    do
        echo -e "./2_openmp_seq $size 1"
        ./3_openmp_seq $size 1
    done
done
