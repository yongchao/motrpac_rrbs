#!/bin/bash
set -euo pipefail
set -x
threads=$1
gdir=$2
shift 2
SID=$(basename $1 _R1.fastq.gz)
paired=0
if(($#==2)); then
    paired=1
    cmd="-1 ../$1 -2 ../$2"
else
    cmd="../$1"
fi    

cd bismark

Ncore=$(($threads/5))
if (($Ncore==0));then
    Ncore=1
fi

#bismark $gdir --ambiguous -un --multicore $Ncore $cmd
bismark $gdir --multicore $Ncore $cmd
ext=_R1_bismark_bt2_
bismark2summary -o ${SID} ${SID}${ext}pe.bam
#mv ${SID}${ext}pe.bam ${SID}.bam
#mv ${SID}${ext}PE_report.txt ${SID}_report.txt

#This was the the one used in nugene, now switced to deduplicate_bismark commmand
#mkdir -p ../tmpdir
#samtools sort -m 2G ${SID}_R1_bismark_bt2_pe.bam  -T ../tmpdir -o ${SID}_sorted.bam -@ $(($threads-1))
#python2 $MOTRPAC_root/nugen/nudup.py -2 --rmdup-only  -o ${SID} -T ../tmpdir ${SID}_sorted.bam
#samtools sort -n -m 2G ${SID}.sorted.dedup.bam -T ../tmpdir -o UMI_${SID}.bam -@ $(($threads-1))
#
#
#bismark_methylation_extractor --multicore $Ncore --comprehensive --bedgraph UMI_${SID}.bam
#        
###remove the temporary files
#rm ${SID}_sorted.bam
#rm ${SID}.sorted.dedup.bam
#
#bSID=UMI_${SID}
#breport=${SID}_R1_bismark_bt2_PE_report.txt
#bismark2report -o ${SID}.html --alignment_report $breport --splitting_report ${bSID}_splitting_report.txt --mbias_report ${bSID}.M-bias.txt

##working with deduplicate_bismark commmand

deduplicate_bismark -p --barcode --bam ${SID}${ext}pe.bam >${SID}_dedup.txt
#debug bismark2summary -o UMI_${SID} ${SID}${ext}pe.deduplicated.bam
#mv ${SID}.deduplicated.bam UMI_${SID}.bam
#mv ${SID}.deduplication_report.txt UMI_${SID}_report.txt
bismark_methylation_extractor --multicore $Ncore --comprehensive --bedgraph ${SID}${ext}pe.deduplicated.bam
bSID=
bismark2report -o ${SID}.html\
	       --dedup_report  ${bSID}_report.txt\
	       --alignment_report ${SID}_report.txt\
	       --splitting_report ${bSID}_splitting_report.txt \
	       --mbias_report ${bSID}.M-bias.txt
