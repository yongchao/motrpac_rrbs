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

#work to do for the lambda to get the bisulfite conversion efficientcy, divided into three groups or not
#
include: "sample_sub.snakefile"
ruleorder: trim>trim_single
        
localrules: bismark_all, pairedness,fastqc_all
rule bismark_all:
    input:        
        expand("bismark/log/OK.{sample}",sample=samples),
        expand("lambda/log/OK.{sample}",sample=samples),
        expand("phix/{sample}.txt",sample=samples),
        expand("chr_info/{sample}.txt",sample=samples),
        "log/OK.pre_align_QC",
        "pairedness"
    output:
        "bismark_qc.csv"
    shell:
        '''
        cd lambda
        set +e #not escaping due to the divison of zero
        bismark2summary
        set -e
        bismark_4strand.sh >bismark_4strand.txt
        
        cd ../bismark
        set +e
	    bismark2summary
        set -e
        bismark_4strand.sh >bismark_4strand.txt

        cd ..
        Rscript --vanilla {root}/bin/qc.R
        '''
rule pairedness:
    output:
        "pairedness"
    run:
        with open(output[0], "w") as out:
            for s in samples:
                out.write(s+"\t"+'%1d'% R2[s]+"\n")
        out.close()

rule trim_single:
    input:        
        trim_input+"{sample}_R1.fastq.gz"
    output:
        "fastq_trim/{sample}_R1.fastq.gz",
        temp("fastq_trim/{sample}_R1_val.fastq.gz")
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
        expand("fastq_trim/{{sample}}_{R}.fastq.gz",R=["R1","R2"]),
        temp(expand("fastq_trim/{{sample}}_{R}_val.fq.gz",R=["R1","R2"]))
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
        state="bismark/log/OK.{sample}",
        bam="bismark/{sample}_R1_bismark_bt2_pe.deduplicated.bam"
    log:
        "bismark/log/{sample}.log"
    threads: 20
    shell:
        '''
        bismark.sh {threads} {gdir} bismark {tmpdir} {input} >&{log}
        echo OK>{output.state}
        '''
#this takes file before MSPI removal
def trim_fastq_info(wildcards):
    sid=wildcards.sample
    return(expand("fastq_trim/{sample}_{R}_val.fq.gz",sample=sid,R=R[sid]))

rule phix:
    input:
        trim_fastq_info
    output:
        "phix/{sample}.txt"
    threads: 6
    shell:
        '''
        gdir_root=$(dirname {gdir})
        gref=$gdir_root/misc_data/phix/phix
        out_tmp=phix/{wildcards.sample}_tmp.txt
        bowtie2.sh -d phix $gref {threads} {input} >& $out_tmp
        mv $out_tmp {output}
        '''
def fastqc_all_input(wildcards):
    files=expand("fastqc/{sample}_fastqc.html",sample=fastqc_samples)
    if fastq_ini=="fastq_raw/":
        files.extend(expand("fastqc_raw/{sample}_fastqc.html",sample=fastqc_samples))
    return(files)    
                
rule fastqc_all:
    input:
        fastqc_all_input
    output:
        "log/OK.fastqc"
    shell:
        '''
        echo OK >{output}      
        '''
rule pre_align_QC:
    input:
        "log/OK.fastqc"
    output:
        "log/OK.pre_align_QC"
    log:
        "log/pre_align_QC.log"
    shell:
        '''
        mkdir -p multiqc
        fastqc_raw=""
        if [[ {fastq_ini} == "fastq_raw/" ]]; then
           fastqc_raw="fastqc_raw fastq_trim"
        fi
        multiqc.sh pre_align fastqc $fastqc_raw >&{log}
        echo OK >{output}
        '''
rule bismark_lambda:
    input:
        fastq_info
    output:
        "lambda/log/OK.{sample}"
    log:
        "lambda/log/{sample}.log"
    threads: 20
    shell:
        '''
        gdir_root=$(dirname {gdir})
        gref=$gdir_root/misc_data/lambda
        bismark.sh {threads} $gref lambda {tmpdir} {input} >&{log}
        echo OK>{output}
        '''
rule chr_info:
    input:
        "bismark/{sample}_R1_bismark_bt2_pe.deduplicated.bam"
    output:
        "chr_info/{sample}.txt"
    shell:
        '''                                                                                                                                                                                                        bam_chrinfo.sh {input}                                                                                                                                                                                     '''


