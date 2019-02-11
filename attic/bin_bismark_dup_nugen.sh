#!/bin/bash
#This is used to test if the bismark gave similar duplcation %, only kept for the historial reference
#The nugen code is fast, but have other possbile problems than bismark
#1. if a read that is mapped to the reverse strand is at an ealier position thatn the read that is mapped to the forward, then this read
#   is ingored
#2. It doens't take the strand the begin and end positions into account  as it only cares about the postions that is mapped to the forward.
set -euo pipefail
set -x
SID=$1
threads=$2
bam=${SID}_R1_bismark_bt2_pe.bam
set +e #head has a problem with this
len=$(samtools view $bam |head -1 |awk '{umi=gensub("^.*:","",1,$1); print length(umi)}')
set -e
samtools sort -m 2G $bam  -T ../tmpdir/ -o ${SID}_sorted.bam -@ $(($threads-1))

#the length maybe wrong as the default is 6, which is no good
python2 $MOTRPAC_root/nugen/nudup.py -2 --rmdup-only  -s $len -l $len -o ${SID} -T ../tmpdir ${SID}_sorted.bam

#rm ${SID}_sorted.bam
#rm ${SID}.sorted.dedup.bam
#rm ${SID}_dup_log.txt
