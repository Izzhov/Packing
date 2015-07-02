#!/bin/bash

dirname=SC-AR-Packings-7-2-15

mkdir ~/Packing/ModSpherocyl/Cluster/$dirname

for j in $(seq 5); do
for i in $(seq 1000); do
    echo cd ~/Packing/ModSpherocyl/Cluster/$dirname '&&' ../ModSpherocyl.out --seed $((1433949739+10000+i)) --ar $((2.0+j*0.1)) 
done
done
