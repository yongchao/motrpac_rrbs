options(stringsAsFactors=FALSE)
#running under the work root folder
getqc<-function(folder){
    x<-read.delim(file.path(folder,"bismark_summary_report.txt"),sep="\t",
                  head=TRUE,row.names=1,check.names=FALSE)
    rownames(x)<-sub("_R1_bismark_bt2_pe.bam$","",rownames(x))
    rownames(x)<-sub("_R1_bismark_bt2.bam$","",rownames(x)) #for single ended read
    y<-read.delim(file.path(folder,"bismark_4strand.txt"),sep="\t",
                  head=TRUE,row.names=1,check.names=FALSE)
    id<-sort(rownames(x))
    x<-cbind(x[id,],y[id,])
    x<-cbind(id=rownames(x),x)
    y<-data.frame(x[,1:2],"%Uniq"=x[,3]/x[,2]*100,
                  "%Unaligned"=x[,4]/x[,2]*100,
                  "%Ambi"=x[,5]/x[,2]*100,
                  "%OT"=x[,16]/rowSums(x[,16:19])*100,
                  "%OB"=x[,17]/rowSums(x[,16:19])*100,
                  "%COT"=x[,18]/rowSums(x[,16:19])*100,
                  "%COB"=x[,19]/rowSums(x[,16:19])*100,
                  "%pos_dup"=x[,7]/x[,3]*100,
                  "%CpG"=x[,10]/(x[,10]+x[,11])*100,
                  "%CHG"=x[,12]/(x[,12]+x[,13])*100,
                  "%CHH"=x[,14]/(x[,14]+x[,15])*100,check.names=FALSE)
    if(dir.exists("fastq_attach")){
        colnames(y)[10]<-"%umi_dup"
    }
    y
}

y<-getqc("bismark")
lambda<-getqc("lambda")
lambda<-lambda[,-c(1,2,4,6:9)]
colnames(lambda)<-paste0("lambda_",colnames(lambda))
#combine these two together
y<-cbind(y,lambda)
readqcinfo<-function(type, name)
{
    read.delim(paste0("multiqc/",type,"_data/multiqc_",name,".txt"),
               sep="\t",head=TRUE,row.names=1,check.names=FALSE,strip.white=TRUE)
}
checknames<-function(x,title)
{
    if(sum(rownames(x)!=samples)==0) return(1)
    #to fix the problem of multqc that the order of S1 is later than S10
    if(sum(sort(rownames(x))!=sort(samples))==0) return(2)
    stop(paste0(title," is wrong"))
}

##read the fastqc info
samples<-y[,1]
rownames(y)<-samples
y<-y[,-1]
fastqc<-readqcinfo("pre_align","fastqc")[,c("Total Sequences","%GC","total_deduplicated_percentage")]
if(2*length(grep("_R2$",rownames(fastqc)))==nrow(fastqc)){
    ##paired, so we do the average of R1 and R2
    id<-(1:(nrow(fastqc)/2))*2
    if(sum(sub("_R2$","_R1",rownames(fastqc)[id])
           !=rownames(fastqc)[id-1])!=0){
        stop("the R1 and R2 files do not match in the fastqc info")
    }
    fastqc<-(fastqc[id-1,]+fastqc[id,])/2
}
rownames(fastqc)<-sub("_R1$","",rownames(fastqc))
#read the trim info and also clean-up fastqc when necessary
TRIM<-dir.exists("fastq_trim")
if(TRIM){
    trim<-readqcinfo("pre_align","cutadapt")
    ##we need to average for the trim_galore
    id<-(1:(nrow(trim)/2))*2
    if(sum(sub("_R2$","_R1",rownames(trim)[id])
           !=rownames(trim)[id-1])!=0){
        stop("the R1 and R2 files do not match in the trim info")
    }
    trim<-(trim[id-1,]+trim[id,])/2
    rownames(trim)<-sub("_R[12]$","",rownames(trim))
    if(checknames(trim,"trim")==2){
        trim<-trim[samples,]
    }
    fastqc_raw<-fastqc[grep("^fastqc_raw \\| ",rownames(fastqc)),]
    rownames(fastqc_raw)<-sub("^fastqc_raw \\| ","",rownames(fastqc_raw))
    if(checknames(fastqc_raw,"fastq_raw")==2){
        fastqc_raw<-fastqc_raw[samples,]
    }
    fastqc<-fastqc[grep("^fastqc \\| ",rownames(fastqc)),]
}
rownames(fastqc)<-sub("^fastqc \\| ","",rownames(fastqc))
if(checknames(fastqc,"fastqc")==2){
    fastqc<-fastqc[samples,]
}
qc<-NULL
if(TRIM) qc<-cbind("reads_raw"=fastqc_raw[,1],"%trimmed"=round(100-as.numeric(y[,1])/fastqc_raw[,"Total Sequences"]*100,dig=2),
                   "%trimmed_bases"=round(trim[,"percent_trimmed"],dig=2))
qc<-cbind(qc,reads=y[,1],"%GC"=round(fastqc[,"%GC"],dig=2),"%dup_sequence"=round(100-fastqc[,"total_deduplicated_percentage"],dig=2))
##read phix
phix<-rep(NA,length(samples))
for(i in 1:length(samples)){
    SID<-samples[i]
    zz<-pipe(paste0("tail -1 ","phix/",SID,".txt |tr -d %"))
    phix[i]<-scan(zz,n=1,quiet=TRUE)
    close(zz)
}
y<-cbind(qc,"%phix"=round(phix,dig=2),round(y[,-1],dig=2))

write.table(y,"bismark_qc.csv",sep=",",col.names=NA,row.name=TRUE)
