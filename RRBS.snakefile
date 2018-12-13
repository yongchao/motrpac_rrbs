#Author: Yongchao Ge

#The usage are in the file RRBS_README.md

#Folder structures at root of working folder
#fastq_raw: raw fastq files, with no adpaters removed, probably a softlink to the fastq files in the output folder of bcl2fastq

#For each sample, we have three possible fastq files
#${sid}_R1.fastq.gz, required for all data. 
#${sid}_R2.fastq.gz, required for pairied end reads
#${sid}_I1.fastq.gz, required for NuGEN with UMI for UMI processsing

#configure genome by using --config genome=hg38_gencode_v29 etc
#Updating the link for hg or hg38 so that we don't have to specify the whole name

include: "sample_sub.snakefile"

localrules: bismark_all
rule bismark_all:
    input:        
        expand("bismark/log/OK.{sample}",sample=samples),
        #expand("qc_insertsize/log/OK.{sample}",sample=samples)
    output:
        "log/OK.bismark"
    shell:
        '''
        cd bismark
	bismark2summary
        echo OK >..{output}
        '''
#rule qc_insertsize:
#    input:
#        "bismark/{sample}_dedup.bam"
#    output:
#        "qc_insertsize/log/OK.{sample}"
#    log:
#        "qc_insertsize/log/{sample}.log"
#    shell:
#        '''
#        qc_insertsize.R {input} qc_insertsize/{wildcards.sample}.pdf >&{log}
#        echo "OK" >{output}
#        '''
rule trim_single:
    input:        
        trim_input+"{sample}_R1.fastq.gz"
    output:
        "fastq_trim/{{sample}}_R1.fastq.gz"
    log:
        "fastq_trim/log/log.{sample}"
    shell:
        '''
        trim_rrbs.sh {input} >&{log}
        '''
rule trim:
    input:        
        expand(trim_input+"{{sample}}_{R}.fastq.gz",R=["R1","R2"])
    output:
        expand("fastq_trim/{{sample}}_{R}.fastq.gz",R=["R1","R2"])
    priority:
        10 #This is preferred than trim_single, with default priority value of zero
    log:
        "fastq_trim/log/log.{sample}"
    shell:
        '''
        trim_rrbs.sh {input} >&{log}
        '''
rule bismark:
    input:
        fastq_info
    output:
        "bismark/log/OK.{sample}"
    log:
        "bismark/log/{sample}.log"
    threads: 20
    shell:
        '''
        bismark.sh {threads} {gdir} {input} >&{log}
        echo OK>{output}
        '''
