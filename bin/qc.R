options(stringsAsFactors=FALSE)
##running under the work root folder
methcap<-scan("methcap.switch")
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
                  "%CTOT"=x[,18]/rowSums(x[,16:19])*100,
                  "%CTOB"=x[,19]/rowSums(x[,16:19])*100,
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
NS<-length(samples)
y<-y[,-1,drop=FALSE]
fastqc<-readqcinfo("pre_align","fastqc")[,c("Total Sequences","%GC","total_deduplicated_percentage")]
PAIRED<-FALSE
if(2*length(grep("_R2$",rownames(fastqc)))==nrow(fastqc)){
    PAIRED<-TRUE
    ##paired, so we do the average of R1 and R2
    id<-(1:(nrow(fastqc)/2))*2
    if(sum(sub("_R2$","_R1",rownames(fastqc)[id])
           !=rownames(fastqc)[id-1])!=0){
        stop("the R1 and R2 files do not match in the fastqc info")
    }
    fastqc<-(fastqc[id-1,,drop=FALSE]+fastqc[id,,drop=FALSE])/2
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
    trim<-(trim[id-1,,drop=FALSE]+trim[id,,drop=FALSE])/2
    rownames(trim)<-sub("_R[12]$","",rownames(trim))
    if(checknames(trim,"trim")==2){
        trim<-trim[samples,]
    }
    fastqc_raw<-fastqc[grep("^fastqc_raw \\| ",rownames(fastqc)),]
    rownames(fastqc_raw)<-sub("^fastqc_raw \\| ","",rownames(fastqc_raw))
    if(checknames(fastqc_raw,"fastq_raw")==2){
        fastqc_raw<-fastqc_raw[samples,,drop=FALSE]
    }
    fastqc<-fastqc[grep("^fastqc \\| ",rownames(fastqc)),,drop=FALSE]
}
rownames(fastqc)<-sub("^fastqc \\| ","",rownames(fastqc))
if(checknames(fastqc,"fastqc")==2){
    fastqc<-fastqc[samples,,drop=FALSE]
}

##read phix and adapter_detected
misc<-matrix(NA,NS,4)
colnames(misc)<-c("phix","adapter_detected","trimmed","no_MSPI")
for(i in 1:NS){
    SID<-samples[i]
    zz<-pipe(paste0("tail -1 ","phix/",SID,".txt |tr -d %"))
    misc[i,1]<-scan(zz,n=1,quiet=TRUE)
    close(zz)

    if(TRIM){
        zz<-pipe(paste0("grep \"^Reads with adapters:\" fastq_trim/log/log.",SID,"|awk -F '[(%]' '{print $2}'"))
        zval<-scan(zz,quiet=TRUE)
        if(length(zval)!= PAIRED+1){
            stop("the fastq_trim log for the contained adapter% is not consistent with the pairedness of the fastq data")
        }
        misc[i,2]<-mean(zval)
        close(zz)

        if (dir.exists("fastq_attach") && methcap==0){##for RRBS
            zz<-pipe(paste0("grep \"Fwd:  D0:\" fastq_trim/log/log.",SID))
            zval<-scan(zz,"",quiet=TRUE)
            close(zz)
            if(length(zval)!=7 && zval[1]!="Fwd:" &&
               length(grep("^other=",zval[6]))!=1 &&
               length(grep("^(total=",zval[7]))!=1)
            {
                stop("The diversity adapter information is incorrect for sample",
                 SID)
            }
            z1<-as.numeric(sub("other=","",zval[6]))
            z2<-as.numeric(sub("\\(total=(.*)\\)","\\1",zval[7],perl=TRUE))
            misc[i,3]<-100-z2/fastqc_raw[i,1]*100
            misc[i,4]<-z1/z2*100
        }else{
            misc[i,3]<-round(100-as.numeric(y[i,1])/fastqc_raw[i,1]*100,dig=3)
        }
    }
}
colnames(misc)<-paste0("%",colnames(misc))
misc<-round(misc,dig=3)


##Read the chr_info.txt
chr_info<-matrix(NA,NS,5)
colnames(chr_info)<-c("chrX","chrY","chrM","chrAuto","contig")

readchr<-function(sid){
    file<-paste0("chr_info/",sid,".txt")
    no.line=scan(pipe(paste0("wc -l ",file)),n=1)
    cat(no.line,"\n")
    if(no.line==0){
        x<-data.frame(c(0,0,0,0))
        rownames(x)<-c("chrX","chrY","chrM","chr1")
    }else{
        x<-read.table(file,
                      row.names=2,strip.white=TRUE)
    }
    xt<-sum(x[,1])
    y<-x[c("chrX","chrY","chrM"),1]/xt*100
    y<-c(y,sum(x[grep("^chr[1-9]",rownames(x)),1])/xt*100)
    y[is.na(y)]<-0
    c(y,100-sum(y))
}
    
for(i in 1:NS){
    chr_info[i,1:5]<-readchr(samples[i])
}
colnames(chr_info)<-paste0("%",colnames(chr_info))
##keep chry with four digits
chr_info[,2]<-round(chr_info[,2],dig=5)
chr_info[,-2]<-round(chr_info[,-2],dig=3)

qc<-NULL
if(TRIM) qc<-cbind("reads_raw"=fastqc_raw[,1],misc[,-1,drop=FALSE],
                   "%trimmed_bases"=round(trim[,"percent_trimmed"],dig=3),
                   "%removed"=round(100-as.numeric(y[,1])/fastqc_raw[,1]*100,dig=3))
qc<-cbind(qc,reads=y[,1],"%GC"=round(fastqc[,"%GC"],dig=3),"%dup_sequence"=round(100-fastqc[,"total_deduplicated_percentage"],dig=3))

y<-cbind(qc,"%phix"=misc[,1],chr_info,round(y[,-1],dig=3))
colnames(y)<-sub("^%","pct_",colnames(y))
colnames(y)<-sub("^lambda_%","lambda_pct_",colnames(y))
if(nrow(y)==1) rownames(y)<-samples
write.table(y,"bismark_qc.csv",sep=",",col.names=NA,row.name=TRUE)
