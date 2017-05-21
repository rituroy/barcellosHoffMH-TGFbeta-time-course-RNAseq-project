## Use Voom with minCnt=10

## ---------------------------------

computerFlag="cluster"
computerFlag=""


if (computerFlag=="cluster") {
    setwd("/home/royr/project/barcellosHoffMH/tgfbTimeCourse")
} else {
    dirSrc="/Users/royr/UCSF/"
    dirSrc2=dirSrc
    setwd(paste(dirSrc2,"barcellosHoffMH/tgfbTimeCourse",sep=""))
}

## ----------------------------------------------

datadir="/Users/royr/Downloads/tmp/"
datadir=""
datadir="results/funcAnno/"
#tbl <- read.table(paste(datadir,"ann_TransgenicMouseVsWildtypeMouse_fromDavid.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)

## -------------------
datadir="results/comparison/"

phen <- read.table(paste("docs/sampleInfo.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
#ann10 <- read.table(paste(datadir,"ann_Homo_sapiens.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
cnt10 <- read.table(paste(datadir,"count_raw_Homo_sapiens.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,row.names="geneId")
#cnt11 <- read.table(paste(datadir,"count_Homo_sapiens.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,row.names="geneId")
phen10 <- read.table(paste(datadir,"sample_Homo_sapiens.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
phen11 <- phen10
stat1_4 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_4hrs_autoDisp.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
stat1_8 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_8hrs_autoDisp.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
stat1_12 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_12hrs_autoDisp.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
stat1_36 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_36hrs_autoDisp.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
#stat1_36 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_36hrs_commonDisp.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
#stat1_36 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_36hrs_tagwiseDisp.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
#stat1_36 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_36hrs_trendedDisp.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)

## -------------------
datadir="results/comparison/countNorm/donorFixedEffect/"

statVF_4_0 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_4hrs_voom_minCnt0.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
statVF_4_2 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_4hrs_voom_minCnt2.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
statVF_4_5 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_4hrs_voom_minCnt5.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
statVF_4_10 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_4hrs_voom_minCnt10.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
statVF_4_20 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_4hrs_voom_minCnt20.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)

statVF_8_0 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_8hrs_voom_minCnt0.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
statVF_8_2 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_8hrs_voom_minCnt2.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
statVF_8_5 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_8hrs_voom_minCnt5.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
statVF_8_10 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_8hrs_voom_minCnt10.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
statVF_8_20 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_8hrs_voom_minCnt20.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)

statVF_12_0 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_12hrs_voom_minCnt0.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
statVF_12_2 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_12hrs_voom_minCnt2.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
statVF_12_5 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_12hrs_voom_minCnt5.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
statVF_12_10 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_12hrs_voom_minCnt10.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
statVF_12_20 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_12hrs_voom_minCnt20.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)

statVF_36_0 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_36hrs_voom_minCnt0.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
statVF_36_2 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_36hrs_voom_minCnt2.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
statVF_36_5 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_36hrs_voom_minCnt5.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
statVF_36_10 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_36hrs_voom_minCnt10.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
statVF_36_20 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_36hrs_voom_minCnt20.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)

## -------------------
datadir="results/comparison/countNorm/donorFixedEffect/"

statVR_4_0 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_4hrs_voom_minCnt0.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
statVR_4_2 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_4hrs_voom_minCnt2.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
statVR_4_5 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_4hrs_voom_minCnt5.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
statVR_4_10 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_4hrs_voom_minCnt10.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
statVR_4_20 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_4hrs_voom_minCnt20.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)

datadir="results/comparison/countNorm/donorRndEffect/"

statVR_8_0 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_8hrs_voom_minCnt0.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
statVR_8_2 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_8hrs_voom_minCnt2.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
statVR_8_5 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_8hrs_voom_minCnt5.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
statVR_8_10 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_8hrs_voom_minCnt10.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
statVR_8_20 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_8hrs_voom_minCnt20.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)

statVR_12_0 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_12hrs_voom_minCnt0.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
statVR_12_2 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_12hrs_voom_minCnt2.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
statVR_12_5 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_12hrs_voom_minCnt5.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
statVR_12_10 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_12hrs_voom_minCnt10.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
statVR_12_20 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_12hrs_voom_minCnt20.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)

statVR_36_0 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_36hrs_voom_minCnt0.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
statVR_36_2 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_36hrs_voom_minCnt2.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
statVR_36_5 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_36hrs_voom_minCnt5.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
statVR_36_10 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_36hrs_voom_minCnt10.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
statVR_36_20 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_36hrs_voom_minCnt20.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)

## -------------------
datadir="results/comparison/countNorm/edgeR/"

statE_4_0 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_4hrs_autoDisp_minCnt0.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
statE_4_2 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_4hrs_autoDisp_minCnt2.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
statE_4_5 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_4hrs_autoDisp_minCnt5.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
statE_4_10 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_4hrs_autoDisp_minCnt10.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
statE_4_20 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_4hrs_autoDisp_minCnt20.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)

statE_8_0 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_8hrs_autoDisp_minCnt0.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
statE_8_2 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_8hrs_autoDisp_minCnt2.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
statE_8_5 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_8hrs_autoDisp_minCnt5.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
statE_8_10 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_8hrs_autoDisp_minCnt10.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
statE_8_20 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_8hrs_autoDisp_minCnt20.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)

statE_12_0 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_12hrs_autoDisp_minCnt0.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
statE_12_2 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_12hrs_autoDisp_minCnt2.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
statE_12_5 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_12hrs_autoDisp_minCnt5.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
statE_12_10 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_12hrs_autoDisp_minCnt10.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
statE_12_20 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_12hrs_autoDisp_minCnt20.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)

statE_36_0 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_36hrs_autoDisp_minCnt0.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
statE_36_2 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_36hrs_autoDisp_minCnt2.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
statE_36_5 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_36hrs_autoDisp_minCnt5.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
statE_36_10 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_36hrs_autoDisp_minCnt10.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
statE_36_20 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_36hrs_autoDisp_minCnt20.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)

## -------------------
stat1_4_0=statVF_4_0
stat1_4_2=statVF_4_2
stat1_4_5=statVF_4_5
stat1_4_10=statVF_4_10
stat1_4_20=statVF_4_20

stat1_8_0=statVF_8_0
stat1_8_2=statVF_8_2
stat1_8_5=statVF_8_5
stat1_8_10=statVF_8_10
stat1_8_20=statVF_8_20

stat1_12_0=statVF_12_0
stat1_12_2=statVF_12_2
stat1_12_5=statVF_12_5
stat1_12_10=statVF_12_10
stat1_12_20=statVF_12_20

stat1_36_0=statVF_36_0
stat1_36_2=statVF_36_2
stat1_36_5=statVF_36_5
stat1_36_10=statVF_36_10
stat1_36_20=statVF_36_20


## -------------------
if (F) {
statE_4=stat1_4
statE_8=stat1_8
statE_12=stat1_12
statE_36=stat1_36

stat1_4$foldChange=2^stat1_4$logFC
stat1_8$foldChange=2^stat1_8$logFC
stat1_12$foldChange=2^stat1_12$logFC
stat1_36$foldChange=2^stat1_36$logFC
#write.table(stat12,file="",col.names=T,row.names=F, sep="\t",quote=F)
#write.table(stat12,file="",col.names=T,row.names=F, sep="\t",quote=F)
#write.table(stat12,file="",col.names=T,row.names=F, sep="\t",quote=F)

#write.table(unique(rownames(cnt10)[rownames(cnt10)%in%tbl[,1]]), file="ensGeneId_TransgenicMouseVsWildtypeMouse_fromDavid.txt",col.names=F,row.names=F, sep="\t",quote=F)
}


cnt_1=t(t(cnt10[,match(phen10$id,colnames(cnt10))])*phen10$norm.factors)
write.table(cbind(geneId=rownames(cnt_1),cnt_1), file="count_norm_Homo_sapiens.txt",col.names=T,row.names=F, sep="\t",quote=F)

if (F) {
    #source(paste(dirSrc,"functions/biomartify.1.1.R",sep=""))
    source(paste(dirSrc,"functions/biomartify.1.2.R",sep=""))
    top=stat1_36[1:10,]
    rownames(top)=top$geneId
    top1 <- biomartify(top, organism = "HomoSapiens")
}

library(biomaRt)

tmpC=rep("",nrow(cnt10))
ann10=data.frame(geneId=rownames(cnt10),geneSym=tmpC,stringsAsFactors=F)
dataset <- "hsapiens_gene_ensembl"
ensembl <- useMart("ensembl", dataset = dataset)
attrs <- c("ensembl_gene_id", "hgnc_symbol", "hgnc_id", "chromosome_name", "start_position", "end_position","strand", "gene_biotype", "description")
attrs <- c("ensembl_gene_id", "hgnc_symbol", "hgnc_id", "entrezgene", "chromosome_name", "start_position", "end_position","strand", "gene_biotype", "description")
res=getBM(attributes=attrs, filters = "ensembl_gene_id", values = rownames(cnt10), mart=ensembl, curl = NULL, checkFilters = TRUE, verbose = FALSE, uniqueRows = TRUE, bmHeader = FALSE)
i=match(rownames(cnt10),res$ensembl_gene_id); i1=which(!is.na(i)); i2=i[i1]; i12=which(is.na(i))
if (length(i12)!=0) {
    tmp=res[1:length(i12),]
    for (k in 1:ncol(tmp)) {
        tmp[,k]=NA
    }
    tmp$ensembl_gene_id=rownames(cnt10)[i12]
    res=rbind(res,tmp)
}
id=unique(res[,1][duplicated(res[,1])])
for (k in 1:length(id)) {
    i=which(res[,1]==id[k])
    if (sum(!duplicated(res$chromosome_name[i]))!=1 | sum(!duplicated(res$start_position[i]))!=1 | sum(!duplicated(res$end_position[i]))!=1) cat(k,"\n")
}
## No replicate gene IDs with different position
i=match(rownames(cnt10),res$ensembl_gene_id); i1=which(!is.na(i)); i2=i[i1]
res=res[i2,]
ann10=res
names(ann10)[match(c("ensembl_gene_id"),names(ann10))]="geneId"
write.table(ann10, file=paste("ann_Homo_sapiens.txt",sep=""),col.names=T,row.names=F, sep="\t",quote=F)
res1=res

stat1_4=stat1_4_10
stat1_8=stat1_8_10
stat1_12=stat1_12_10
stat1_36=stat1_36_10
for (compFlag in c("TGFbetaVuntreated_4hrs","TGFbetaVuntreated_8hrs","TGFbetaVuntreated_12hrs","TGFbetaVuntreated_36hrs")) {
    switch(compFlag,
    "TGFbetaVuntreated_4hrs"={
        compName="4hr: TGFbeta vs. untreated"
        compFlag2=compFlag
        colGeneId="geneId"
        stat_1=stat1_4
        ann_1=ann10
    },
    "TGFbetaVuntreated_8hrs"={
        compName="8hr: TGFbeta vs. untreated"
        compFlag2=compFlag
        colGeneId="geneId"
        stat_1=stat1_8
        ann_1=ann10
    },
    "TGFbetaVuntreated_12hrs"={
        compName="12hr: TGFbeta vs. untreated"
        compFlag2=compFlag
        colGeneId="geneId"
        stat_1=stat1_12
        ann_1=ann10
    },
    "TGFbetaVuntreated_36hrs"={
        compName="36hr: TGFbeta vs. untreated"
        compFlag2=compFlag
        colGeneId="geneId"
        stat_1=stat1_36
        ann_1=ann10
    }
    )
    i=1:10
    i=1:nrow(stat_1)
    stat_1=stat_1[i,]
    res=res1
    i=match(stat_1$geneId,res$ensembl_gene_id); i1=which(!is.na(i)); i2=i[i1]; i12=which(is.na(i))
    if (length(i12)!=0) {
        tmp=res[1:length(i12),]
        for (k in 1:ncol(tmp)) {
            tmp[,k]=NA
        }
        tmp$ensembl_gene_id=stat_1$geneId[i12]
        res=rbind(res,tmp)
    }
    i=match(stat_1$geneId,res$ensembl_gene_id); i1=which(!is.na(i)); i2=i[i1]
    names(stat_1)[match("FDR",names(stat_1))]="Qvalue"
    stat_1=cbind(stat_1,res[i2,which(names(res)!="ensembl_gene_id")])
    write.table(stat_1, file=paste("stat_",compFlag2,".txt",sep=""),col.names=T,row.names=F, sep="\t",quote=F)
}


## -------------------
stat1=statVR_36_0; stat2=statE_36_0
stat1=statVR_4_0; stat2=statE_4_0
stat1=statVR_36_20; stat2=statVF_36_20
stat1=statVR_36_20; stat2=statE_36_20
stat1=statVF_36_20; stat2=statE_36_20
i=match(stat1$geneId,stat2$geneId); i1=which(!is.na(i)); i2=i[i1]
lim=range(c(stat1$logFC[i1],stat2$logFC[i2]))
plot(stat1$logFC[i1],stat2$logFC[i2],xlim=lim,ylim=lim)
abline(c(0,1),lty="dotted")

for (timeFlag in c(4,8,12,36)) {
    png(paste("tmp_",timeFlag,".png",sep=""))
    par(mfrow=c(2,2))
    switch(as.character(timeFlag),
    "4"={
        statE=statE_4_10
        statVF=statVF_4_10
        statVR=statVR_4_10
    },
    "8"={
        statE=statE_8_10
        statVF=statVF_8_10
        statVR=statVR_8_10
    },
    "12"={
        statE=statE_12_10
        statVF=statVF_12_10
        statVR=statVR_12_10
    },
    "36"={
        statE=statE_36_10
        statVF=statVF_36_10
        statVR=statVR_36_10
    }
    )
    for (compId in c("eVvf","eVvr","vfVvr")) {
        switch(compId,
        "eVvf"={
         stat1=statE; stat2=statVF
        },
        "eVvf"={
            stat1=statE; stat2=statVR
        },
        "vfVvr"={
            stat1=statVF; stat2=statVR
        }
        )
        header=paste(timeFlag," hr",sep="")
        xlab=strsplit(compId,"V")[[1]][1]
        ylab=strsplit(compId,"V")[[1]][2]
        i=match(stat1$geneId,stat2$geneId); i1=which(!is.na(i)); i2=i[i1]
        lim=range(c(stat1$logFC[i1],stat2$logFC[i2]))
        plot(stat1$logFC[i1],stat2$logFC[i2],xlim=lim,ylim=lim,main=header,xlab=xlab,ylab=ylab)
        abline(c(0,1),lty="dotted")
    }
    dev.off()
}


## -------------------
stat1=statVF_36_10; stat2=statVF_36_20
stat1=statVF_4_10; stat2=statVF_4_20
i1=1:nrow(stat1); i2=1:nrow(stat2)
i=match(stat1$geneId[i1], stat2$geneId[i2]); i1=i1[which(!is.na(i))]; i2=i2[i[i1]]
par(mfrow=c(2,2))
plot(stat1$logFC[i1], stat2$logFC[i2])
plot(stat1$PValue[i1], stat2$PValue[i2])


## ----------------------------------------------
stat_1=statVF_8_10
stat_1=statVF_12_10
stat_1=statVF_36_10
stat_1=statVF_4_10
ann_1=ann10[match(stat_1$geneId,ann10$geneId),]
i=grep("tgfb",tolower(ann_1$hgnc_symbol))
i=i[which(stat_1$PValue[i]<0.05)]
cbind(ann_1[i,c("hgnc_symbol","chromosome_name","start_position")],stat_1[i,][order(stat_1$PValue[i]),])

## -------------------
library(edgeR)


## ----------------------------------------------
## Compare voom & edgeR results

datadir="results/comparison/"
colList=c("red","blue","green")
for (compFlag in c("TGFbetaVuntreated_4hrs","TGFbetaVuntreated_8hrs","TGFbetaVuntreated_12hrs","TGFbetaVuntreated_36hrs")) {
    load(paste(datadir,"dge_Homo_sapiens_",compFlag,".RData",sep=""))
	switch(compFlag,
	"TGFbetaVuntreated_4hrs"={
        stat_1=statVF_36_10; stat_2=statE_36_10; cnt_1=t(t(dgeT$counts)*dgeT$samples$norm.factors); phen_1=cbind(phen[match(rownames(dgeT$samples),phen$id),],dgeT$samples[,c("lib.size","norm.factors")]); ann_1=ann10[match(rownames(dgeT$counts),ann10$geneId),]
        colGeneId="geneId"
	},
    "TGFbetaVuntreated_8hrs"={
        stat_1=statVF_36_10; stat_2=statE_36_10; cnt_1=t(t(dgeT$counts)*dgeT$samples$norm.factors); phen_1=cbind(phen[match(rownames(dgeT$samples),phen$id),],dgeT$samples[,c("lib.size","norm.factors")]); ann_1=ann10[match(rownames(dgeT$counts),ann10$geneId),]
        colGeneId="geneId"
    },
    "TGFbetaVuntreated_12hrs"={
        stat_1=statVF_36_10; stat_2=statE_36_10; cnt_1=t(t(dgeT$counts)*dgeT$samples$norm.factors); phen_1=cbind(phen[match(rownames(dgeT$samples),phen$id),],dgeT$samples[,c("lib.size","norm.factors")]); ann_1=ann10[match(rownames(dgeT$counts),ann10$geneId),]
        colGeneId="geneId"
    },
    "TGFbetaVuntreated_36hrs"={
        stat_1=statVF_36_10; stat_2=statE_36_10; cnt_1=t(t(dgeT$counts)*dgeT$samples$norm.factors); phen_1=cbind(phen[match(rownames(dgeT$samples),phen$id),],dgeT$samples[,c("lib.size","norm.factors")]); ann_1=ann10[match(rownames(dgeT$counts),ann10$geneId),]
        colGeneId="geneId"
    }
	)
	i=match(stat_1[,colGeneId],stat_2[,colGeneId]); i1=which(!is.na(i)); i2=i[i1]
	dat1=stat_1[i1,]; dat2=stat_2[i2,]; ann_2=ann_1[match(dat1[,colGeneId],ann_1[,colGeneId]),]
	for (k in 1:ncol(dat1)) {
		if (any(is.na(dat1[,k])!=is.na(dat2[,k])) | any(dat1[,k]!=dat2[,k],na.rm=T)) {
			cat("\n\n=========",k,names(dat1)[k],"=============\n")
			if (is.numeric(dat1[,k])) {
				png(paste("plot_",names(dat1)[k],"_voomVsEdgeR_for",compFlag,".png",sep=""))
				lim=range(dat1[,k],dat2[,k],na.rm=T)
				plot(dat1[,k],dat2[,k],xlim=lim,ylim=lim,main=compFlag,xlab=paste("Voom: ",names(dat1)[k],sep=""),ylab=paste("edgeR: ",names(dat1)[k],sep=""))
				abline(c(0,1),col="red",lty="dotted")
				dev.off()
				if (names(dat1)[k]%in%c("PValue","FDR")) {
					pThres=0.05
					x=table(voom=dat1[,k]<pThres,edgeR=dat2[,k]<pThres,dnn=c(paste("voom: ",names(dat1)[k],"<",pThres,sep=""),paste("edgeR: ",names(dat1)[k],"<",pThres,sep="")))
					print(x)
					print(summary(dat1[dat2[,k]<pThres,k]))
					print(summary(dat2[dat1[,k]<pThres,k]))
				}
			}
		}
	}

	## ------------------------------

#	dat0=t(t(cnt_1[match(stat_1[,colGeneId][i1],rownames(cnt_1)),])*phen_1$norm.factors)
	dat0=cnt_1[match(stat_1[,colGeneId][i1],rownames(cnt_1)),]
	dat0=dat0+min(c(dat0[dat0!=0]))*0.1
	dat1=log2(dat0)
	#dat1[!is.finite(dat1)]=NA
	#dat1=dat0

    geneName=paste(ann_2$geneId,", ",ann_2$geneSym,sep="")
    ttl=c("tgfb","untr")

    pThres=0.05
    pThres1=0.01; pThres2=0.2
    pThres1=0.05; pThres2=0.05
    ii=order(stat_1$PValue[i1])
    ii=ii[which(stat_1$FDR[i1][ii]<pThres1)]
    set.seed(5345)
    ii=ii[sample(1:length(ii),replace=F)][1:32]
    for (iiii in 1:2) {
        png(paste("plot_count_signifVoom_for",compFlag,"_",iiii,".png",sep=""),width=4*240, height=2*240)
        par(mfcol=c(2,8))
        par(mar=c(5, 4, 4, 2) + 0.1)
        par(mar=c(6, 4, 4, 2) + 0.1)
        for (iii in seq(iiii,length(ii),by=2)) {
            i=ii[iii]
            nm=paste(geneName[i],"\nPV ",signif(stat_1$PValue[i1][i],2),sep="")
            nm=paste(geneName[i],"\nPV ",signif(stat_1$PValue[i1][i],2),"\nlog2FC ",round(stat_1$logFC[i1][i],2)," logCPM ",round(stat_1$logCPM[i1][i],2),sep="")
            stripchart(dat1[i,]~phen_1$treat,xlim=c(.5,2.5),group.names=ttl,sub=nm,ylab="log2(normalized count)",vertical=T,type="n")
            points(as.integer(as.factor(phen_1$treat)),dat1[i,],col=colList[phen_1$donor])
        }
        par(mfcol=c(1,1))
        title(main="Significant genes")
        dev.off()
    }
    
    pThres=0.05
    pThres1=0.01; pThres2=0.2
    pThres1=0.5; pThres2=0.5
    ii=order(stat_1$PValue[i1])
    ii=ii[which(stat_1$FDR[i1][ii]>=pThres1)]
    set.seed(5345)
    ii=ii[sample(1:length(ii),replace=F)][1:32]
    for (iiii in 1:2) {
        png(paste("plot_count_notsignifVoom_for",compFlag,"_",iiii,".png",sep=""),width=4*240, height=2*240)
        par(mfcol=c(2,8))
        par(mar=c(6, 4, 4, 2) + 0.1)
        for (iii in seq(iiii,length(ii),by=2)) {
            i=ii[iii]
            nm=paste(geneName[i],"\nPV ",signif(stat_1$PValue[i1][i],2),sep="")
            nm=paste(geneName[i],"\nPV ",signif(stat_1$PValue[i1][i],2),"\nlog2FC ",round(stat_1$logFC[i1][i],2)," logCPM ",round(stat_1$logCPM[i1][i],2),sep="")
            stripchart(dat1[i,]~phen_1$treat,xlim=c(.5,2.5),group.names=ttl,sub=nm,ylab="log2(normalized count)",vertical=T,type="n")
            points(as.integer(as.factor(phen_1$treat)),dat1[i,],col=colList[phen_1$donor])
        }
        par(mfcol=c(1,1))
        title(main="Not significant genes")
        dev.off()
    }

	pThres=0.05
	pThres1=0.01; pThres2=0.2
	pThres1=0.05; pThres2=0.05
	ii=order(stat_1$PValue[i1])
	ii=ii[which(stat_1$FDR[i1][ii]<pThres1 & stat_2$PValue[i2][ii]>=pThres2)][1:6]
	#ii=ii[which(stat_1$FDR[i1][ii]<pThres1 & stat_2$FDR[i2][ii]>=pThres2)][1:6]
	#ii=ii[order(stat_2$PValue[i2][ii],decreasing=T)][1:6]
	png(paste("plot_count_signifVoomNotsignifEdgeR_for",compFlag,".png",sep=""))
	par(mfcol=c(2,3))
	for (iii in 1:length(ii)) {
		i=ii[iii]
		stripchart(dat1[i,]~phen_1$treat,xlim=c(.5,2.5),sub=paste(geneName[i],"\nPV ",signif(stat_1$PValue[i1][i],2)," (voom), ",signif(stat_2$PValue[i2][i],2)," (edgeR)",sep=""),ylab="log2(normalized count)",vertical=T,type="n")
        points(as.integer(as.factor(phen_1$treat)),dat1[i,],col=colList[phen_1$donor])
	}
	par(mfcol=c(1,1))
	title(main="Signif for voom, not signif for edgeR")
	dev.off()

	pThres1=0.05; pThres2=0.05
	ii=order(stat_2$PValue[i2])
	ii=ii[which(stat_2$FDR[i2][ii]<pThres2 & stat_1$PValue[i1][ii]>=pThres1)][1:6]
	png(paste("plot_count_notsignifVoomSignifEdgeR_for",compFlag,".png",sep=""))
	par(mfcol=c(2,3))
	for (iii in 1:length(ii)) {
		i=ii[iii]
		stripchart(dat1[i,]~phen_1$treat,xlim=c(.5,2.5),sub=paste(geneName[i],"\nPV ",signif(stat_1$PValue[i1][i],2)," (voom), ",signif(stat_2$PValue[i2][i],2)," (edgeR)",sep=""),ylab="log2(normalized count)",vertical=T,type="n")
        points(as.integer(as.factor(phen_1$treat)),dat1[i,],col=colList[phen_1$donor])
	}
	par(mfcol=c(1,1))
	title(main="Not signif for voom, signif for edgeR")
	dev.off()

	pThres1=0.05; pThres2=0.05
	ii=order(stat_1$PValue[i1])
	ii=ii[which(stat_1$FDR[i1][ii]<pThres1 & stat_2$FDR[i2][ii]<pThres2)][1:6]
	png(paste("plot_count_signifVoomEdgeR_for",compFlag,".png",sep=""))
	par(mfcol=c(2,3))
	for (iii in 1:length(ii)) {
		i=ii[iii]
		stripchart(dat1[i,]~phen_1$treat,xlim=c(.5,2.5),sub=paste(geneName[i],"\nPV ",signif(stat_1$PValue[i1][i],2)," (voom), ",signif(stat_2$PValue[i2][i],2)," (edgeR)",sep=""),ylab="log2(normalized count)",vertical=T,type="n")
        points(as.integer(as.factor(phen_1$treat)),dat1[i,],col=colList[phen_1$donor])
	}
	par(mfcol=c(1,1))
	title(main="Signif for both voom and edgeR")
	dev.off()

	pThres1=0.5; pThres2=0.5
	#ii=order(stat_1$PValue[i1],decreasing=T)
	ii=1:length(i1)
	ii=ii[which(stat_1$FDR[i1][ii]>=pThres1 & stat_2$FDR[i2][ii]>=pThres2)][1:6]
	png(paste("plot_count_notsignifVoomEdgeR_for",compFlag,".png",sep=""))
	par(mfcol=c(2,3))
	for (iii in 1:length(ii)) {
		i=ii[iii]
		stripchart(dat1[i,]~phen_1$treat,xlim=c(.5,2.5),sub=paste(geneName[i],"\nPV ",signif(stat_1$PValue[i1][i],2)," (voom), ",signif(stat_2$PValue[i2][i],2)," (edgeR)",sep=""),ylab="log2(normalized count)",vertical=T,type="n")
        points(as.integer(as.factor(phen_1$treat)),dat1[i,],col=colList[phen_1$donor])
	}
	par(mfcol=c(1,1))
	title(main="Not signif for both voom and edgeR")
	dev.off()
}

## ----------------------------------------------
## DE genes

colIdPV="pv"; colNamePV="PV"; pThres=10^-6
colIdPV="qv"; colNamePV="QV"; pThres=0.05

table(stat1_4$PValue<pThres)
table(stat1_8$PValue<pThres)
table(stat1_12$PValue<pThres)
table(stat1_36$PValue<pThres)

onePlotFlag=F
onePlotFlag=T

minCntFlag=c(10)
compFlag3=""

minCntFlag=c(0,2,5,10,20)
compFlag3="_edgerVsVoom"
compFlag3="_edgeR"
compFlag3="_voom"

out=NULL
if (compFlag3=="") {
    if (onePlotFlag) {
        if (compFlag3=="") {
            png(paste("dePlots.png",sep=""),width=4*240,height=2*240)
            par(mfcol=c(2,4))
        } else {
            png(paste("dePlots.png",sep=""),width=6*240,height=2*240)
            par(mfcol=c(2,6))
        }
    }
}
for (compFlag in c("TGFbetaVuntreated_4hrs","TGFbetaVuntreated_8hrs","TGFbetaVuntreated_12hrs","TGFbetaVuntreated_36hrs")) {
    if (compFlag3=="") compList=compFlag
    if (compFlag3=="_edgerVsVoom") compList=c("edgeR",paste("voom_minCnt",minCntFlag,sep=""))
    if (compFlag3=="_voom") compList=paste("voom_minCnt",minCntFlag,sep="")
    if (compFlag3=="_edgeR") compList=paste("edgeR_minCnt",minCntFlag,sep="")
    if (compFlag3!="") {
        if (onePlotFlag) {
            if (compFlag3=="") {
                png(paste("dePlots",compFlag3,"_",compFlag,".png",sep=""),width=4*240,height=2*240)
                par(mfcol=c(2,4))
            } else {
                png(paste("dePlots",compFlag3,"_",compFlag,".png",sep=""),width=length(compList)*240,height=2*240)
                par(mfcol=c(2,length(compList)))
            }
        }
    }
    for (compThis in compList) {
        switch(compThis,
        "edgeR"={
            stat1_4=statE_4
            stat1_8=statE_8
            stat1_12=statE_12
            stat1_36=statE_36
        },
        "edgeR_minCnt0"={
            stat1_4=statE_4_0
            stat1_8=statE_8_0
            stat1_12=statE_12_0
            stat1_36=statE_36_0
        },
        "edgeR_minCnt2"={
            stat1_4=statE_4_2
            stat1_8=statE_8_2
            stat1_12=statE_12_2
            stat1_36=statE_36_2
        },
        "edgeR_minCnt5"={
            stat1_4=statE_4_5
            stat1_8=statE_8_5
            stat1_12=statE_12_5
            stat1_36=statE_36_5
        },
        "edgeR_minCnt10"={
            stat1_4=statE_4_10
            stat1_8=statE_8_10
            stat1_12=statE_12_10
            stat1_36=statE_36_10
        },
        "edgeR_minCnt20"={
            stat1_4=statE_4_20
            stat1_8=statE_8_20
            stat1_12=statE_12_20
            stat1_36=statE_36_20
        },
        "voom_minCnt0"={
            stat1_4=stat1_4_0
            stat1_8=stat1_8_0
            stat1_12=stat1_12_0
            stat1_36=stat1_36_0
        },
        "voom_minCnt2"={
            stat1_4=stat1_4_2
            stat1_8=stat1_8_2
            stat1_12=stat1_12_2
            stat1_36=stat1_36_2
        },
        "voom_minCnt5"={
            stat1_4=stat1_4_5
            stat1_8=stat1_8_5
            stat1_12=stat1_12_5
            stat1_36=stat1_36_5
        },
        "voom_minCnt10"={
            stat1_4=stat1_4_10
            stat1_8=stat1_8_10
            stat1_12=stat1_12_10
            stat1_36=stat1_36_10
        },
        "voom_minCnt20"={
            stat1_4=stat1_4_20
            stat1_8=stat1_8_20
            stat1_12=stat1_12_20
            stat1_36=stat1_36_20
        }
        )
        switch(compFlag,
               "TGFbetaVuntreated_4hrs"={
                   compName="4hr: TGFbeta vs. untreated"
                   compFlag2=compFlag
                   colGeneId="geneId"
                   stat_1=stat1_4
                   ann_1=ann10
               },
               "TGFbetaVuntreated_8hrs"={
                   compName="8hr: TGFbeta vs. untreated"
                   compFlag2=compFlag
                   colGeneId="geneId"
                   stat_1=stat1_8
                   ann_1=ann10
               },
               "TGFbetaVuntreated_12hrs"={
                   compName="12hr: TGFbeta vs. untreated"
                   compFlag2=compFlag
                   colGeneId="geneId"
                   stat_1=stat1_12
                   ann_1=ann10
               },
               "TGFbetaVuntreated_36hrs"={
                   compName="36hr: TGFbeta vs. untreated"
                   compFlag2=compFlag
                   colGeneId="geneId"
                   stat_1=stat1_36
                   ann_1=ann10
               }
               )
        if (compFlag3!="") compName=paste(compThis,", ",compName,sep="")
        i1=match(stat_1[,colGeneId],ann_1[,colGeneId])
        ann_1=ann_1[i1,]
        cat("\n\n",compName,"\n",sep="")
        x=table(stat_1$logFC>0,stat_1$PValue<pThres)
        rownames(x)=c("Down","Up")[match(rownames(x),c("FALSE","TRUE"))]
        colnames(x)=paste(colNamePV,c(">=","<"),pThres,sep="")[match(colnames(x),c("FALSE","TRUE"))]
        print(x)
        
        if (F) {
            lim=c(0,1)
            png(paste("tmp_",fName,".png",sep=""))
            plot(stat$pvBeta[i],stat$qvBeta[i],xlim=lim,ylim=lim,xlab="P-value",ylab="Q-value",pch=19,cex.axis=1.5,cex.lab=1.5,main=compName)
            abline(0,1)
            dev.off()
            
            png(paste("qqplot_",fName,".png",sep=""))
            pvs <- sort(na.exclude(stat$pvBeta[i]))
            qqplot(-log10(runif(length(pvs),0,1)),-log10(pvs),xlab="Expected -log10(p-values) by random",ylab="Observed -log10(p-values)",pch=19,cex.axis=1.5,cex.lab=1.5,main=compName)
            abline(0,1)
            dev.off()
        }
        
        if (!onePlotFlag) png(paste("histogram_",compFlag,".png",sep=""))
        hist(stat_1$PValue,xlab="P-value",pch=19,cex.axis=1.5,cex.lab=1.5,main=compName)
        if (!onePlotFlag) dev.off()
        
        if (!onePlotFlag) png(paste("volcanoPlot_",compFlag,"_",tolower(colNamePV),pThres,".png",sep=""))
        if (colIdPV=="pv") i=which(stat_1$PValue<pThres) else i=which(stat_1$FDR<pThres)
        plot(stat_1$logFC,-log10(stat_1$PValue),xlab="Log2 fold change",ylab="-Log10(p-value)",main=paste(compName,"\nNo. with ",tolower(colNamePV)," < ",pThres,": ",length(i),sep=""),pch=19,cex=.8,cex.axis=1.5,cex.lab=1.5)
        if (length(i)!=0) points(stat_1$logFC[i],-log10(stat_1$PValue[i]),col="red",pch=19,cex=.8)
        if (!onePlotFlag) dev.off()

    #    out2=stat_1[match(ann_1[,colGeneId],stat_1[,colGeneId]),c("foldChange","logFC","logCPM","PValue","FDR")]
        out2=stat_1[match(ann_1[,colGeneId],stat_1[,colGeneId]),c("logFC","PValue")]
        names(out2)=paste(names(out2),"_",compFlag2,sep="")
        out2=cbind(ann_1[,which(!names(ann_1)%in%c("gene_id"))],out2)
        for (k in 1:ncol(out2)) {
            if (is.character(out2[,k])) {
                out2[,k]=gsub(", ...","",out2[,k])
            }
        }
        write.table(out2, file=paste("stat_",compFlag2,".txt",sep=""),col.names=T,row.names=F, sep="\t",quote=F)
        
        write.table(unique(ann_1$geneId), file=paste("ensGeneId_",compFlag,".txt",sep=""),col.names=F,row.names=F, sep="\t",quote=F)
        i=which(stat_1$PValue<pThres)
        if (length(i)!=0) write.table(unique(ann_1$geneId[i]), file=paste("ensGeneId_pv",pThres,"_",compFlag,".txt",sep=""),col.names=F,row.names=F, sep="\t",quote=F)
        i=order(stat_1$PValue)[1:200]
        write.table(unique(ann_1$geneId[i]), file=paste("ensGeneId_top200_",compFlag,".txt",sep=""),col.names=F,row.names=F, sep="\t",quote=F)
    }
    if (onePlotFlag & compFlag3!="") dev.off()
}
if (onePlotFlag & compFlag3=="") dev.off()

## Check log2 fold change direction
j1=which(phen12$type=="untreated")
j2=which(phen12$type=="TGFbeta")
log2Expr=cbind(wt=apply(log2(cnt12[i1,j1]),1,mean,na.rm=T),tr=apply(log2(cnt12[i1,j2]),1,mean,na.rm=T))
i=which(stat_1$PValue<pThres)
cbind(stat_1$logFC[i],log2Expr[i,])[which(stat_1$logFC[i]>0),][1:10,]
cbind(stat_1$logFC[i],log2Expr[i,])[which(stat_1$logFC[i]<0),][1:10,]

###########################################################
###########################################################


## ----------------------------------------------
library(clusterProfiler)

outFormat="png"
outFormat="pdf"

outFormat2="pdf"

compList="TGFbetaVuntreated_8hrs"
compList=c("TGFbetaVuntreated_4hrs","TGFbetaVuntreated_8hrs","TGFbetaVuntreated_12hrs","TGFbetaVuntreated_36hrs")

colIdEst="logFC"; colIdPV=c("PValue","FDR"); pThres=0.05
colInfo=data.frame(key=c("PValue","FDR"),value=c("pv","qv"),stringsAsFactors=F)
#for (pThres in c(0.05)) {
for (pThres in c(NA)) {
    for (compFlag in compList) {
        switch(compFlag,
        "TGFbetaVuntreated_4hrs"={
            compName="4hr: TGFbeta vs. untreated"
            compFlag2=compFlag
            colGeneId="geneId"
            stat_1=stat1_4
            ann_1=ann10
        },
        "TGFbetaVuntreated_8hrs"={
            compName="8hr: TGFbeta vs. untreated"
            compFlag2=compFlag
            colGeneId="geneId"
            stat_1=stat1_8
            ann_1=ann10
        },
        "TGFbetaVuntreated_12hrs"={
            compName="12hr: TGFbeta vs. untreated"
            compFlag2=compFlag
            colGeneId="geneId"
            stat_1=stat1_12
            ann_1=ann10
        },
        "TGFbetaVuntreated_36hrs"={
            compName="36hr: TGFbeta vs. untreated"
            compFlag2=compFlag
            colGeneId="geneId"
            stat_1=stat1_36
            ann_1=ann10
        }
        )
        cat("\n\n================== ",compFlag,", ",pThres," ==================\n",sep="")
        for (pThres1 in c(0.05,0.01,0.001,0.0005)) {
            x=table(stat_1$FDR<pThres1)
            names(x)=paste("QV",c(">=","<"),pThres1,sep="")[match(names(x),c("FALSE","TRUE"))]
            print(x)
        }
        
        
        iA2=match(stat_1$geneId,ann_1$geneId)
        
        idType(annoDb = "org.Hs.eg.db")
        out=bitr(unique(ann_1$hgnc_symbol[!is.na(ann_1$hgnc_symbol)]), fromType="SYMBOL", toType=c("ENTREZID", "ENSEMBL"), annoDb="org.Hs.eg.db")
        table(is.na(as.integer(out$ENTREZID))) ## ALL FALSE
        ann_1$entrezId=""
        i=match(ann_1$hgnc_symbol,out$SYMBOL); i1=which(!is.na(i)); i2=i[i1]
        ann_1$entrezId[i1]=out$ENTREZID[i2]
        table(ann_1$entrezId=="")
        
        ann2=ann_1[iA2,]
        if (is.na(pThres)) {
            fName2=paste("_",compFlag,"_top200genes",sep="")
            prId=order(stat_1[,colIdPV[2]])[1:200]
        } else {
            fName2=paste("_",compFlag,"_",colInfo$value[which(colInfo$key==colIdPV[2])],pThres,sep="")
            prId=which(stat_1[,colIdPV[2]]<pThres)
        }
        
        #for (enrichFlag in c("go","kegg","do","david")) {
        #for (enrichFlag in c("go","kegg")) {
        for (enrichFlag in c("go")) {
            ontList=""
            switch(enrichFlag,
            "go"={
                ontList=c("BP","CC","MF")
                pThres2=0.05; qThres2=0.05
            },
            "kegg"={
                pThres2=0.4; qThres2=0.4
                pThres2=0.05; qThres2=0.05
            },
            "do"={
                pThres2=0.05; qThres2=0.05
            },
            "david"={
            }
            )
            pThres2=0.005; qThres2=0.05
            pThres2=0.05; qThres2=0.05
            
            for (ontThis in ontList) {
                cat("\n\n================== ",compFlag,", ",pThres,", ",enrichFlag,", ",ontThis," ==================\n",sep="")
                i0=which(!duplicated(ann2$entrezId) & ann2$entrezId!="")
                i=prId[which(!duplicated(ann2$entrezId[prId]) & ann2$entrezId[prId]!="")]
                res=NULL
                switch(enrichFlag,
                "go"={
                    cat("------------- groupGO\n\n",sep="")
                    fName3=paste("groupGO",fName2,ifelse(ontThis=="","","_"),ontThis,sep="")
                    res=NULL
                    #res=groupGO(gene=ann2$entrezId[prId],organism="human",ont=ontThis,level=3,readable=TRUE)
                    res=groupGO(gene=ann2$entrezId[i],organism="human",ont=ontThis,level=3,readable=TRUE)
                    #print(head(summary(res)))
                    if (outFormat=="png") {
                        png(paste(fName3,".png",sep=""))
                    } else if (outFormat=="pdf") {
                        pdf(paste(fName3,".pdf",sep=""))
                    }
                    #plot(1:5,1:5)
                    barplot(res, drop=TRUE, showCategory=12)
                    dev.off()
                    if (F) {
                        png(paste("tmp",fName2,ifelse(ontThis=="","","_"),ontThis,".png",sep=""))
                        plot(1:5,1:5)
                        dev.off()
                    }
                    res=NULL
                    res=try(enrichGO(gene=ann2$entrezId[i],universe=ann2$entrezId[i0],organism="human",ont=ontThis,pAdjustMethod="BH",pvalueCutoff=pThres2,qvalueCutoff=qThres2,readable=TRUE))
                },
                "kegg"={
                    res=enrichKEGG(gene=ann2$entrezId[i],universe=ann2$entrezId[i0],pAdjustMethod="BH",pvalueCutoff=pThres2,qvalueCutoff=qThres2,readable=TRUE,use_internal_data=TRUE)
                },
                "do"={
                    res=enrichDO(gene=ann2$entrezId[i],ont="DO",universe=ann2$entrezId[i0],pAdjustMethod="BH",pvalueCutoff=pThres2,qvalueCutoff=qThres2,minGSSize=5,readable=FALSE)
                },
                "david"={
                    res=enrichDAVID(gene=ann2$entrezId[i],idType="ENTREZ_GENE_ID",listType="Gene",annotation="KEGG_PATHWAY",david.user="clusterProfiler@hku.hk")
                }
                )
                head(summary(res))
                if (class(res) %in% c("try-error","NULL")) {
                    print(paste(ontThis,": Error",sep=""))
                } else {
                    if (nrow(res@result)==0) {
                        print(paste(ontThis,": Nothing significant",sep=""))
                    } else {
                        fName3=paste(enrichFlag,"Anno",ifelse(ontThis=="","","_"),ontThis,"_pv",pThres2,"_qv",qThres2,fName2,sep="")
                        write.table(res@result, file=paste(fName3,".txt",sep=""),col.names=T,row.names=F, sep="\t",quote=F)
                        cat("PV<",signif(max(res@result$pvalue,na.rm=T),2),", QV<",signif(max(res@result$qvalue,na.rm=T),2),"\n")
                        fName3=paste(enrichFlag,"EnrichMap",ifelse(ontThis=="","","_"),ontThis,"_pv",pThres2,"_qv",qThres2,fName2,sep="")
                        if (outFormat=="png") {
                            png(paste(fName3,".png",sep=""))
                        } else if (outFormat=="pdf") {
                            pdf(paste(fName3,".pdf",sep=""), pointsize=4)
                        }
                        try(enrichMap(res,vertex.label.font=.4))
                        dev.off()
                        fName3=paste(enrichFlag,"CgetPlot",ifelse(ontThis=="","","_"),ontThis,"_pv",pThres2,"_qv",qThres2,fName2,sep="")
                        if (outFormat=="png") {
                            #png(paste(fName3,".png",sep=""),width=480,height=480,pointsize=12)
                            png(paste(fName3,".png",sep=""),width=1.5*480,height=1.5*480,pointsize=12)
                        } else if (outFormat=="pdf") {
                            pdf(paste(fName3,".pdf",sep=""))
                        }
                        try(cnetplot(res, categorySize="pvalue", foldChange=2^stat_1$coef[i0]))
                        dev.off()
                        if (enrichFlag=="go") {
                            fName3=paste(enrichFlag,"PlotGraph",ifelse(ontThis=="","","_"),ontThis,"_pv",pThres2,"_qv",qThres2,fName2,sep="")
                            if (outFormat2=="png") {
                                #png(paste(fName3,".png",sep=""),width=480,height=480,pointsize=12)
                                #png(paste(fName3,".png",sep=""),width=480,height=480,pointsize=24)
                                png(paste(fName3,".png",sep=""),width=5*480,height=5*480,pointsize=24)
                            } else if (outFormat2=="pdf") {
                                pdf(paste(fName3,".pdf",sep=""))
                            }
                            try(plotGOgraph(res))
                            dev.off()
                        }
                    }
                }
            }
        }
    }
}

###########################################################
###########################################################

## Run heatmap.R

###########################################################
###########################################################
