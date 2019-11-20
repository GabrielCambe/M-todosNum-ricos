#!/bin/bash

for g in FLOPS_DP; do
    for n in $1; do
	likwid-perfctr -O -C 3 -g $g -m ./matmult -n $n > debug$g.log
	raw_data=$(cat debug$g.log);
	for region in multMatPtrVet; do 
	    metrica=$(echo $raw_data | awk -F'[:,;]' '{ for(i=1; i<=NF; i++){ if(($i~"TABLE")&&($(i+1)~"Region multMatPtrVet")&&($(i+2)~"Group 1 Metric")&&($(i+3)~"FLOPS_DP")){ for(j=i+1; j<=NF; j++){ if($j~"TABLE"){ exit; } else if($j~"DP MFLOP\/s"){ print $(j+1); }}}}}')
	    echo $n $metrica >> $g$region.dat
	done

	for region in multMatRowVet; do 
	    metrica=$(echo $raw_data | awk -F'[:,;]' '{ for(i=1; i<=NF; i++){ if(($i~"TABLE")&&($(i+1)~"Region multMatRowVet")&&($(i+2)~"Group 1 Metric")&&($(i+3)~"FLOPS_DP")){ for(j=i+1; j<=NF; j++){ if($j~"TABLE"){ exit; } else if($j~"DP MFLOP\/s"){ print $(j+1); }}}}}')
	    echo $n $metrica >> $g$region.dat
	done

	for region in multMatColVet; do 
	    metrica=$(echo $raw_data | awk -F'[:,;]' '{ for(i=1; i<=NF; i++){ if(($i~"TABLE")&&($(i+1)~"Region multMatColVet")&&($(i+2)~"Group 1 Metric")&&($(i+3)~"FLOPS_DP")){ for(j=i+1; j<=NF; j++){ if($j~"TABLE"){ exit; } else if($j~"DP MFLOP\/s"){ print $(j+1); }}}}}')
	    echo $n $metrica >> $g$region.dat
	done
    done
done
