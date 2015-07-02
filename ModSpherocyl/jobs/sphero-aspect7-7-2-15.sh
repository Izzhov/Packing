#!/bin/bash

dirname=SC-AR-Packings-7-2-15

mkdir ~/Packing/ModSpherocyl/Cluster/$dirname

for j in $(seq 8); do
for i in $(seq 1000); do
    echo cd ~/Packing/ModSpherocyl/Cluster/$dirname '&&' ../ModSpherocyl.out --seed $((1433949739+22000+i)) --ar $((3.2+j*0.1)) 
done
done
