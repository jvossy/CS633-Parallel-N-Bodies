#!/bin/bash
FILENAME=bench.dat 
for i in 10 20 50 100 200 500 1000 2000 5000 10000
do
    ./v1 $i >> $FILENAME
    echo >> $FILENAME
done