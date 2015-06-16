#!/bin/bash

dirname=6-2D-Packings-6-16-15

mkdir ~/Packing/ModSph2DSteep/Cluster/$dirname

for i in $(seq 100000); do
    echo cd ~/Packing/ModSph2DSteep/Cluster/$dirname '&&' ../Sphere2DSteepest.out --seed $((1433949739+i)) --filename $i.txt
done
