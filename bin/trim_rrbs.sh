#!/bin/bash -x
set -eu -o pipefail

index_adapter=AGATCGGAAGAGC
univ_adapter=AAATCAAAAAAAC
#This program include two steps, one is trim_galore
#and the othr is nugen way of removing diversity

R1=$(basename $1 .fastq.gz)
fd=$(dirname $1)
if [[ $fd == fastq_attach ]]; then
    umi=1
else
    umi=0
fi
       
R2=""
opt2=""
paired=0
if (( $# == 2 )); then
    R2=$(basename $2 .fastq.gz)
    opt2="--paired -a2 $univ_adapter"
    paired=1
fi

trim_galore -a $index_adapter $opt2 -o fastq_trim "$@"
cd fastq_trim
#first need to do the clean-up

#Trim_galore outputs different output depending on the paired or non-paired
if (( $paired == 0)); then
    mv ${R1}_trimmed.fq.gz ${R1}_val_1.fq.gz
fi

if [[ $(dirname $1) ==  fastq_attach ]]; then
    #UMI operations
    #need extra step for removing the diversity
    opt2=""
    if(($paired==1));then
	opt2="-2 ${R2}_val_2.fq.gz"
    fi
    python2 $MOTRPAC_root/nugen/trimRRBSdiversityAdaptCustomers.py -1 ${R1}_val_1.fq.gz $opt2
    mv ${R1}_val_1.fq_trimmed.fq.gz $R1.fastq.gz
    rm ${R1}_val_1.fq.gz
    if(($paired==1));then
	mv ${R2}_val_2.fq_trimmed.fq.gz $R2.fastq.gz
	rm ${R2}_val_2.fq.gz
    fi
else
    mv ${R1}_val_1.fq.gz ${R1}.fastq.gz
    if (( $paired == 1)); then
	mv ${R2}_val_2.fq.gz ${R2}.fastq.gz
    fi

fi
