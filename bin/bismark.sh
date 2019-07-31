#!/bin/bash
set -euo pipefail
set -x
threads=$1
gdir=$2
odir=$3
tmpdir_root=$4
shift 4
SID=$(basename $1 _R1.fastq.gz)
if(($#==2)); then
    ##paired
    cmd="-1 ../$1 -2 ../$2"
    optp="-p"
    align_report=${SID}_R1_bismark_bt2_PE_report.txt
    bam=${SID}_R1_bismark_bt2_pe.bam
else
    #single
    cmd="../$1"
    optp="-s"
    align_report=${SID}_R1_bismark_bt2_SE_report.txt
    bam=${SID}_R1_bismark_bt2.bam
fi    

cd $odir

Ncore=$(($threads/4))
if (($Ncore==0));then
    Ncore=1
fi
date

tmpdir=$(mktemp -d -p $tmpdir_root bismark.${SID}.XXX)
bismark $gdir --multicore $Ncore --temp_dir $tmpdir $cmd

#sort the bam files according to the best alignment score and then by the readname
#and write the sorted bam file back to the original file
sam=$tmpdir/$SID.sam
samtools view -@ $threads -H $bam >$sam
samtools view -@ $threads $bam |sort -k5,5nr -k1,1 -s -S30G -T $tmpdir >>$sam
samtools view -@ $threads -b $sam -o $bam

rm -rf $tmpdir

#find out how many have been matched
nreads=$(samtools view -c $bam)

if(($nreads==0)); then
    #no-deduplicates
    exit 0
fi

if [[ -e ../fastq_attach ]]; then
    #UMI
    date
    bismark_bam_UMI_format.sh $bam
    date
    deduplicate_bismark $optp --barcode --bam $bam >${SID}_dedup.txt
else
    #Just use the position to deduplicate, good for MethCAP data
    deduplicate_bismark $optp --bam $bam >${SID}_dedup.txt
fi
bam=$(basename $bam .bam).deduplicated.bam
date
bismark_methylation_extractor --multicore $Ncore --comprehensive --bedgraph $bam
bismark2report -o ${SID}.html\
	       -a $align_report

