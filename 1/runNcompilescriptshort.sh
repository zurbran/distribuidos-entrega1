#!/bin/bash

cd 2
cd openmp
make
make sequential

for size in 512 1024 2048
do
    for threads in 2 4
    do
        for runs in 1 2 3 4
        do
            echo -e "./2_openmp $size $threads"
            ./2_openmp $size $threads
        done
    done
done

for size in 512 1024 2048
do
    for runs in 1 2 3 4
    do
        echo -e "./2_openmp_seq $size 1"
        ./2_openmp_seq $size 1
    done
done
