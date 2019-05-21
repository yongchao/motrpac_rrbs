#!/bin/bash
set -euo pipefail
bam=$1
SID=$(basename $bam _R1_bismark_bt2_pe.deduplicated.bam)
samtools view $bam|cut -f 3|sort |uniq -c > chr_info/${SID}.txt
