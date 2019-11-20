#!/bin/bash

echo "performance" > /sys/devices/system/cpu/cpufreq/policy3/scaling_governor

rm -f *.dat

#for g in L3 L2CACHE FLOPS_DP FLOPS_AVX; do
for g in FLOPS_AVX; do
#for n in 32 50 64 100 128 200 256 300 400 512 1000 1024 2000 2048 3000 4000 4096 5000 10000; do
    for i in 32 50; do
	./medir-perf$g.sh $i
    done
    ./$g\plot.gplt
done

echo "powersave" > /sys/devices/system/cpu/cpufreq/policy3/scaling_governor
