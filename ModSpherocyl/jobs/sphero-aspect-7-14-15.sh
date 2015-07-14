#!/bin/bash

dirname=SC-AR-Packings-7-14-15

mkdir ~/Packing/ModSpherocyl/Cluster/$dirname

for k in $(seq 3); do
for j in $(seq 9); do
for i in $(seq 1000); do
    echo cd ~/Packing/ModSpherocyl/Cluster/$dirname '&&' ../ModSpherocyl.out --seed $((1435762516+(k-1)*10000+(j-1)*1000+(i-1))) --filename $((1435762516+(k-1)*10000+(j-1)*1000+(i-1))).txt --ar $k.$j
done
done
for i in $(seq 1000); do
    echo cd ~/Packing/ModSpherocyl/Cluster/$dirname '&&' ../ModSpherocyl.out --seed $((1435762516+k*10000+(i-1))) --filename $((1435762516+(k-1)*10000+(j-1)*1000+(i-1))).txt --ar $((k+1)).0
done
done
