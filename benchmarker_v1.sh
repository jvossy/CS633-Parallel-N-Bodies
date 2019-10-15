#!/bin/bash
FILENAME=data/sequentialv1_bench.csv
rm $FILENAME
touch $FILENAME

gcc v1.c -fopenmp -Ofast -o v1

for i in 10 20 50 100 200 500 1000 2000 5000 10000
do
    simple_sum=0
    for iter in {1..10}
    do
        simple_sum=$(echo "${simple_sum}+$(./v1 $i )"    | bc -l)
    done
    simple_average=$(echo "${simple_sum} / 10.0"    | bc -l)
    echo "${i},${simple_average}" >> $FILENAME
done
