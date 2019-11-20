#!/bin/bash

echo "performance" > /sys/devices/system/cpu/cpufreq/policy3/scaling_governor

rm *.dat

for g in L2CACHE; do
#    for n in 32 50 64 100 128 200 256 300 400 512 1000 1024 2000 2048 3000 4000 4096 5000 10000; do
    for n in 32 50 64 100 128; do
	likwid-perfctr -O -C 3 -g $g -m ./matmult -n $n > debug.log
	raw_data=$(cat debug.log);
	for region in multMatPtrVet; do 
	    metrica=$(echo $raw_data | awk -F'[:,;]' '{ for(i=1; i<=NF; i++){ if(($i~"TABLE")&&($(i+1)~"Region multMatPtrVet")&&($(i+2)~"Group 1 Metric")&&($(i+3)~"L2CACHE")){ for(j=i+1; j<=NF; j++){ if($j~"TABLE"){ exit; } else if($j~/L2 miss ratio/){ print $(j+1); }}}}}')
	    echo $n $metrica >> $g$region.dat
	done

	for region in multMatRowVet; do 
	    metrica=$(echo $raw_data | awk -F'[:,;]' '{ for(i=1; i<=NF; i++){ if(($i~"TABLE")&&($(i+1)~"Region multMatRowVet")&&($(i+2)~"Group 1 Metric")&&($(i+3)~"L2CACHE")){ for(j=i+1; j<=NF; j++){ if($j~"TABLE"){ exit; } else if($j~/L2 miss ratio/){ print $(j+1); }}}}}')
	    echo $n $metrica >> $g$region.dat
	done

	for region in multMatColVet; do 
	    metrica=$(echo $raw_data | awk -F'[:,;]' '{ for(i=1; i<=NF; i++){ if(($i~"TABLE")&&($(i+1)~"Region multMatColVet")&&($(i+2)~"Group 1 Metric")&&($(i+3)~"L2CACHE")){ for(j=i+1; j<=NF; j++){ if($j~"TABLE"){ exit; } else if($j~/L2 miss ratio/){ print $(j+1); }}}}}')
	    echo $n $metrica >> $g$region.dat
	done
    done
done

echo "powersave" > /sys/devices/system/cpu/cpufreq/policy3/scaling_governor
