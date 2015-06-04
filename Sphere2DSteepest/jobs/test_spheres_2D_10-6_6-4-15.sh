#!/bin/bash

dirname=testsphere2D_10-6_6-4-15

mkdir ~/Packing/Sphere2DSteepest/Cluster/$dirname

for i in $(seq 100); do
    echo cd ~/Packing/Sphere2DSteepest/Cluster/$dirname '&&' ../Sphere2DSteepest.out --seed $((1433355833+i)) --filename $i.txt
done
