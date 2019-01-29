#!/bin/bash
set -euo pipefail
#running under the bismark sub folder 
sample=(`cut -f 1 ../pairedness`)
pairedness=(`cut -f 2 ../pairedness`)
#ext can be _R1_bismark_bt2_pe.deduplicated_splitting_report.txt or  _R1_bismark_bt2.deduplicated_splitting_report.txt
printf "samples\treads\ttotal C\tmeth CpG\tmeth CHG\tmeth CHH\tunmeth CPG\tunmeth CHG\tunmeth CHH\n"
for i in ${!sample[@]}; do
    s=${sample[$i]}
    ext=""
    if ((${pairedness[$i]}==1)); then
	ext="_pe"
    fi
    ext=_R1_bismark_bt2$ext.deduplicated_splitting_report.txt
    
    printf $s"\t"
    grep -P "Processed|Total" ${s}$ext | grep ":" | cut -f 2 -d":"|tr -d " \t"| tr "\n" "\t"|sed 's/\t$//'
    printf "\n"
done
