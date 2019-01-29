# snakemake implementation of MoTrPAC RRBS/MethCAP pipeline
**Contact:** Yongchao Ge (yongchao.ge@mssm.edu)
* [The MoTrPAC RRBS MOP (web view version 2.0)](https://docs.google.com/document/d/e/2PACX-1vTxMwrq4Q3b5GUfPtZF2krpK_ah0yW--TyeOAFrEVi_FvgIPhkCKPytRQ8QmZe5WF1KKjah0pftU9A_/pub)
* [The MoTrPAC MethCAP MOP (web view version 2.0)](https://docs.google.com/document/d/e/2PACX-1vT_qPrhekYh8VMDVy3ACGYapnTol6aUmekR6-zh_10RR0jLXiUkfse9Y6KyTuMS2KDpOnoeEPM8mbVC/pub)
* [rnq-seq snakemake implementation](https://github.com/yongchao/motrpac_rnaseq)

Please note this pipeline was primarily written for RRBS pipeline, but it also works for MethCAP pipeline with different setting of the fastq files. For the RRBS pipeline, it required three fastq files (`${SID}_R1.fastq.gz`, `${SID}_I1.fastq.gz` and `${SID}_R2.fastq.gz`) in the `fastq_raw` folder, while for MethCAP pipeline, it required two fastq files (`${SID}_R1.fastq.gz`, and `${SID}_R2.fastq.gz`) in the fastq_raw folder. The pipeline will testing the absence of `_I1` file to switch to the MethCAP pipeline (see details in Section B below).

# A. External softwares installation and bash environmental setup
## A.1 install extra software required for RRBS data
This snakemake implementation assumes that you already have set-up the [rna-seq snakemake pipeline](https://github.com/yongchao/motrpac_rnaseq).
In addition to the software used in rna-seq data, we need to install additional software by issuing the command below. The details are seen in [conda\_install\_extra.sh](bin/conda_install_extra.sh)
```bash
conda install \
      seqtk=1.3 \
      trim-galore=0.5.0
```
## A.2 Git clone the code for the motrpac_rrbs
As this snakemake implementation uses the code that is part of rna-seq snakemake pipeline, the folder position to where to put the code for morpac-rrbs is important.
The command `git clone https://github.com/yongchao/motrpac_rrbs` needs to be issued under the parent folder of `motrpac\_rnaseq` folder that contains the git code of rna-seq snakemake pipeline. The resulting two folders `motrpac\_rrbs` (contains the git code for RRBS snakemake pipeline) and `motrpac_rnaseq` (contains the git code for the RNA-seq snakemake pipeline) should have the same parent folder. This ensure all of the softlinks in the `motrpac\_rnaseq` folder can function properly.

## A.3 Bash environment setup
Following the instructions from [RNA-seq README file section A.2](https://github.com/yongchao/motrpac_rnaseq/blob/master/README.md#a2-bash-environments-setup) to export the environmental variables `PATH`, `MOTRPAC_root`,`MOTRPAC_conda`,`MOTRPAC_refdata`

Please be advised at the `MOTRPAC_root` needs to be pointed to the git code base of the rrbs, not rna-seq, so you should run the command 
`export $(motrapac_rrbs/bin/load_motrpac.sh-c $conda -r $refdata)`  not `export $(motrapac_rnseq/bin/load_motrpac.sh-c $conda -r $refdata)`, i.e., the file load_motrpac.sh should come from RRBS git code base rather than the rna-seq git code base.

## A.4 Download the genome source and build the refdata
This is already done by the rnq-seq pipeline, we don't need to do anything more

## B.2 Run the snakemake program
* In a work folder, a subfolder `fastq_raw` contains the fastq files of all samples `${SID}_R1.fastq.gz`, `${SID}_I1.fastq.gz` and `${SID}_R2.fastq.gz`.
* If the `${SID}_I1.fastq.gz` is missing, we will assume the data is for MethCAP data. (Please note for the folder `fasstq_raw`, either all `_I1` files are missing or all `_I1` files exist)
* Make sure the `MOTRPAC_root`, `PATH` and other environmental variables have been setup correctly according to section A.2
* Run the command locally to debug possible problems below for the human genome  
  `snakemake -s $MOTRPAC_root/RRBS.snakefile`
* If the data is for rat samples, run the command below for the rat genome `rn6_ensembl_r95`  
  `snakemake -s $MOTRPAC_root/RRBS.snakefile --config genome=rn6_ensembl_r95`
* If the snakemake is running OK locally, then submit the snakemake jobs to the cluster. This is only necessary for large jobs. This script was written for Sinai LSF jobs submission system. Other cluster job submission system may need to write their own script.  
  `Snakemake.lsf -- -s $MOTRPAC_root/RRBS.snakefile --config genome=rn6_ensembl_r95`
  
