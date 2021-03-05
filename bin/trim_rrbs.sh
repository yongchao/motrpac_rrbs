#!/bin/bash -x
set -eu -o pipefail
methcap=$1
shift

if [[ $methcap == 0 ]]; then 
    index_adapter=AGATCGGAAGAGC
    index2_adapter=AAATCAAAAAAAC
else
    index_adapter=AGATCGGAAGAGC
    index2_adapter=$index_adapter
fi
#This program include two steps, one is trim_galore
#and the othr is nugen way of removing diversity

R1=$(basename $1 .fastq.gz)
fd=$(dirname $1)
if [[ $fd == fastq_attach ]]; then
    umi=1
else
    umi=0
    index2_adapter=$index_adapter #note that with no umi, we are using the MethCAP data and with MethCAP primer
fi
       
R2=""
opt2=""
paired=0
if (( $# == 2 )); then
    R2=$(basename $2 .fastq.gz)
    opt2="--paired -a2 $index2_adapter"
    paired=1
fi

trim_galore -a $index_adapter $opt2 -o fastq_trim "$@"
cd fastq_trim
#first need to do the clean-up

#Trim_galore outputs different output depending on the paired or non-paired
if (( $paired == 0)); then
    mv ${R1}_trimmed.fq.gz ${R1}_val_1.fq.gz
fi

if [[ $(dirname $1) ==  fastq_attach && $methcap == 0 ]]; then
    #UMI operations
    #need extra step for removing the diversity
    opt2=""
    if(($paired==1));then
	opt2="-2 ${R2}_val_2.fq.gz"
    fi
    python2 $MOTRPAC_root/nugen/trimRRBSdiversityAdaptCustomers.py -1 ${R1}_val_1.fq.gz $opt2
    mv ${R1}_val_1.fq_trimmed.fq.gz $R1.fastq.gz
    mv ${R1}_val_1.fq.gz ${R1}_val.fq.gz
    if(($paired==1));then
	mv ${R2}_val_2.fq_trimmed.fq.gz $R2.fastq.gz
	mv ${R2}_val_2.fq.gz ${R2}_val.fq.gz
    fi
else
    mv ${R1}_val_1.fq.gz ${R1}.fastq.gz
    if (( $paired == 1)); then
	mv ${R2}_val_2.fq.gz ${R2}.fastq.gz
    fi

fi
