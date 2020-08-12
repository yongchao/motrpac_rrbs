#!/bin/bash
set -euo pipefail
bam=$1
if [[ $bam == *_pe.deduplicated.bam ]]; then
    SID=$(basename $bam _R1_bismark_bt2_pe.deduplicated.bam)
else
    SID=$(basename $bam _R1_bismark_bt2.deduplicated.bam)
fi

samtools view $bam|cut -f 3|sort |uniq -c > chr_info/${SID}.txt
