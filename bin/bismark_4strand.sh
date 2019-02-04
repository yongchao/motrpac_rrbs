#!/bin/bash
set -euo pipefail
#This is under bismark subfolder
sample=(`cut -f 1 ../pairedness`)
pairedness=(`cut -f 2 ../pairedness`)
printf "samples\tOT\tOB\tCTOT\tCTOB\n"
for i in ${!sample[@]}; do
    s=${sample[$i]}
    ext="_SE"
    if ((${pairedness[$i]}==1)); then
	ext="_PE"
    fi
    ext=_R1_bismark_bt2${ext}_report.txt
    printf $s"\t"
    grep "((converted) top strand)" -A 3 ${s}$ext |sort | cut -f 2 |tr "\n" "\t"|sed 's/\t$//'
    #The sort is necessary due to the different ordering between single and paied ends for the
    #four strands
    printf "\n"
done
