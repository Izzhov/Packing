#!/bin/bash

dirname=SC-AR-Packings-7-14-15

mkdir ~/Packing/ModSpherocyl/Cluster/$dirname

b=$(( ($1-1)*1500 ))
c=$((b+1500))

for k in $(seq 3); do
for j in $(seq 10); do
for i in $(seq 1000); do
    a=$(( (k-1)*10000+(j-1)*1000+(i-1) ))
    b=$(( ($1-1)*1500 ))
    c=$((b+1500))
    if ((a >= b)); then
    if ((a < c)); then
    echo cd ~/Packing/ModSpherocyl/Cluster/$dirname '&&' ../ModSpherocyl.out --seed $((1435762516+a)) --filename $((1435762516+a)).txt --ar $k.$j
    fi
    fi
done
done
done
