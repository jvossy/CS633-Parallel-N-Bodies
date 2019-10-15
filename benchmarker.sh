#!/bin/bash
FILENAME=bench.dat
rm $FILENAME
gcc v1par.c -fopenmp -Ofast -o v1par
for i in 10 20 50 100 200 500 1000 2000 5000 10000
do
    ./v1par $i >> $FILENAME
    echo >> $FILENAME
done