FILENAME=data/sequentialv3_bench.dat
rm $FILENAME
touch $FILENAME

gcc v3.c -fopenmp -Ofast -o v3

for i in 10 20 50 100 200 500 1000 2000 5000 10000
do
    simple_sum=0
    for iter in {1..10}
    do
        simple_sum=$(echo "${simple_sum}+$(./v3 $i )"    | bc -l)
    done
    simple_average=$(echo "${simple_sum} / 10.0"    | bc -l)
    echo "${i}, ${simple_average}" >> $FILENAME
done

