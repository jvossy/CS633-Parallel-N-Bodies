#!/bin/bash
SIMPLE_FILENAME=data/parallelv1_bench.csv
rm $SIMPLE_FILENAME
touch $SIMPLE_FILENAME

gcc v1par.c -fopenmp -Ofast -o v1par

for i in 200 500 1000 2000 5000 10000 20000
do
    for num_threads in 1 4 8 16 18 20
    do
        sum=0
        export OMP_NUM_THREADS=$num_threads
        for iter in {1..10}
        do
            sum=$(echo "${sum}+$(./v1par $i )"    | bc -l)
        done
        average=$(echo "${sum} / 10.0"    | bc -l)
        echo "${i}, ${num_threads}, ${average}" >> $SIMPLE_FILENAME
    done
done
