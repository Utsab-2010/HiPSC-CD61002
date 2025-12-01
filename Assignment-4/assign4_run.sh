#!/bin/bash

mpicc assgn4.c -o assgn4 -lm

OUT="results_assgn4.csv"
echo "N,P,Time,Iterations,Residual" > $OUT

OMEGA=1.9
TOL=1e-6
MAX_ITER=10000

N_LIST=(80 120 200)
P_LIST=(1 2 8 16 32)

for n in "${N_LIST[@]}"; do
    echo "Running N=$n with Omega=$OMEGA"
    for p in "${P_LIST[@]}"; do
        result=$(mpirun  -np $p ./assgn4 $n $OMEGA $TOL $MAX_ITER)
        echo "$result" >> $OUT
        echo "  P=$p | $result"
    done
done