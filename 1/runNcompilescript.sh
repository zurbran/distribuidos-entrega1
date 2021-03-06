#!/bin/bash

cd 1
cd openmp
make
make sequential

for size in 512 1024 2048
do
    for threads in 2 4
    do
        for runs in 1 2 3 4
        do
            echo -e "./1_openmp $size $threads"
            ./1_openmp $size $threads
        done
    done
done

for size in 512 1024 2048
do
    for runs in 1 2 3 4
    do
        echo -e "./1_openmp_seq $size 1"
        ./1_openmp_seq $size 1
    done
done

cd ..
cd pthreads
make

for size in 512 1024 2048
do
    for threads in 2 4
    do
        for runs in 1 2 3 4
        do
            echo -e "./1_pthreads $size $threads"
            ./1_pthreads $size $threads
        done
    done
done

cd ..
cd ..
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

cd ..
cd pthreads
make

for size in 512 1024 2048
do
    for threads in 2 4
    do
        for runs in 1 2 3 4
        do
            echo -e "./2_pthreads $size $threads"
            ./2_pthreads $size $threads
        done
    done
done

cd ..
cd ..
cd 3
cd openmp

make
make sequential

for size in "160000000" "320000000" "640000000"
do
    for threads in 2 4
    do
        for runs in 1 2 3 4
        do
            echo -e "./3_openmp $size $threads"
            ./3_openmp $size $threads
        done
    done
done

for size in "160000000" "320000000" "640000000"
do
    for runs in 1 2 3 4
    do
        echo -e "./3_openmp_seq $size 1"
        ./3_openmp_seq $size 1
    done
done
