#!/bin/bash
SIMPLE_FILENAME=data/simple_bench.dat
NEWTON_FILENAME=data/newton_bench.dat
rm $SIMPLE_FILENAME
rm $NEWTON_FILENAME
touch $SIMPLE_FILENAME
touch $NEWTON_FILENAME

gcc v1par.c -fopenmp -Ofast -o v1par
gcc v3npar.c -fopenmp -Ofast -o v3npar

for i in 10 20 50 100 200 500 1000 2000 5000 10000
do
    simple_sum=0
    newton_sum=0
    for iter in {1..10}
    do
        simple_sum=$(echo "${simple_sum}+$(./v1par $i )"    | bc -l)
        newton_sum=$(echo "${newton_sum}+$(./v3npar $i)"    | bc -l)
    done
    simple_average=$(echo "${simple_sum} / 10.0"    | bc -l)
    newton_average=$(echo "${newton_sum} / 10.0"    | bc -l)
    echo "${i}     ${simple_average}" >> $SIMPLE_FILENAME
    echo "${i}     ${newton_average}" >> $NEWTON_FILENAME
done
