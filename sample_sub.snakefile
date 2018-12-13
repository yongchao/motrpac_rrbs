#This file is shared between RRBS.snakefile and RNA-seq.snakefile
root=os.environ['MOTRPAC_root']
if "genome" in config:
    gdir=os.environ['MOTRPAC_refdata']+"/"+config["genome"]
else:
    gdir=os.environ['MOTRPAC_refdata']+"/hg38_gencode_v29"

#If softlins fastq_raw to fastq folder, then no trim is happening
    
#do not allow to go to the sub directory
wildcard_constraints:
    sample="[^/]+"
    
#If fastq_raw exists, then trim and all folders are based on fastq_trim
#Otherwise lookinf for files in fastq folder
import os
if os.path.isdir("fastq"):
    fastq_ini="fastq/"
    fastq="fastq/"
elif os.path.isdir("fastq_raw"):
    fastq_ini="fastq_raw/"
    fastq="fastq_trim/"
else:
    print("The sub folders fastq_raw or fastq do not exit, exit\n")
    sys.exit(1)


samples,=glob_wildcards(fastq_ini+"{sample,[^/]+}_R1.fastq.gz")
fastqc_samples,=glob_wildcards(fastq_ini+"{sample,[^/]+_R[12]}.fastq.gz")

#construct R2 and I1 info
R2={}
R={}
I=0    
for s in samples:
    if os.path.isfile(fastq_ini+s+"_R2.fastq.gz"):
        R2[s]=1
        R[s]=["R1","R2"]
    else:
        R2[s]=0
        R[s]=["R1"]
        
    if os.path.isfile(fastq_ini+s+"_I1.fastq.gz"):
        I+=1

samples_all=[s for s in samples] #this is for summarization, it may include UMI
if I>0:
    if I!=len(samples):
        print("Not all samples have UMI index files, exit\n")
        print(samples)
        sys.exit(1)
    else:
        trim_input="fastq_attach/"
        for s in samples:
            samples_all.extend(["UMI_"+s])
            R2["UMI_"+s]=R2[s] #not necessary, just to avoid keyerror
            R["UMI_"+s]=R[s]
        
else:
    trim_input="fastq_raw/"

    
#All samples should have I or not
#We could mix smaples with single ends and paired ends
        
if(len(samples)==0):
    print("There are no fastq files in fastq_raw or fastq folders, exit\n")
    sys.exit(1)
        
def R2_info(wildcards):
    return(R2[wildcards.sample])
    
def fastq_info(wildcards,prefix=fastq):
    sid=wildcards.sample
    return(expand(prefix+"{sample}_{R}.fastq.gz",sample=sid,R=R[sid]))

ruleorder: trim > trim_single
rule samples_all:
    output:
        "samples_all",
        "samples"
    run:
        with open(output[0], "w") as out:
            for s in samples_all:
                out.write(s+"\n")
        out.close()        
        with open(output[1], "w") as out:
            for s in samples:
                out.write(s+"\n")
        out.close()
                
rule fastqc: #on the final fastq data before star
    input:
        fastq+"{sample}.fastq.gz"
    output:
        "fastqc/{sample}_fastqc.html"
    log:
        "fastqc/log/{sample}.log"
    shell:
        '''
        fastqc.sh {input} fastqc >& {log}
        '''
        
rule fastqc_raw: #on the raw data
    input:
        "fastq_raw/{sample}.fastq.gz"
    output:
        "fastqc_raw/{sample}_fastqc.html"
    log:
        "fastqc_raw/log/{sample}.log"
    shell:
        '''
        fastqc.sh {input} fastqc_raw >& {log}
        '''
        
rule UMI_attach:
    input:
        "fastq_raw/{sample}.fastq.gz"
    output:
        "fastq_attach/{sample}.fastq.gz"
    log:
        "fastq_attach/log/{sample}.log"
    shell:
        '''
        SID={wildcards.sample}
        I=${{SID%R[12]}}I1  #remove R[12] and then add back to I1
        zcat {input} | UMI_attach.awk -v Ifq=fastq_raw/$I.fastq.gz|gzip -c>{output} 
        '''
