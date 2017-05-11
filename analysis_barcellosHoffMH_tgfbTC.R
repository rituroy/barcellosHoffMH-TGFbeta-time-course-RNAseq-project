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

## -------------------
datadir="results/final/tables/"

statVF_4_0 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_4hrs_voom_minCnt0.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
statVF_4_2 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_4hrs_voom_minCnt2.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
statVF_4_5 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_4hrs_voom_minCnt5.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
statVF_4_10 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_4hrs.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
statVF_4_20 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_4hrs_voom_minCnt20.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)

statVF_8_0 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_8hrs_voom_minCnt0.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
statVF_8_2 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_8hrs_voom_minCnt2.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
statVF_8_5 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_8hrs_voom_minCnt5.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
statVF_8_10 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_8hrs.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
statVF_8_20 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_8hrs_voom_minCnt20.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)

statVF_12_0 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_12hrs_voom_minCnt0.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
statVF_12_2 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_12hrs_voom_minCnt2.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
statVF_12_5 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_12hrs_voom_minCnt5.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
statVF_12_10 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_12hrs.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
statVF_12_20 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_12hrs_voom_minCnt20.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)

statVF_36_0 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_36hrs_voom_minCnt0.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
statVF_36_2 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_36hrs_voom_minCnt2.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
statVF_36_5 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_36hrs_voom_minCnt5.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
statVF_36_10 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_36hrs.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
statVF_36_20 <- read.table(paste(datadir,"stat_TGFbetaVuntreated_36hrs_voom_minCnt20.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)

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


stat1_4=statVF_4_10
stat1_8=statVF_8_10
stat1_12=statVF_12_10
stat1_36=statVF_36_10


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
source("biomartify.R")
source("funcs1.R")
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

minCntFlag=c(0,2,5,10,20)
compFlag3="_edgerVsVoom"
compFlag3="_edgeR"
compFlag3="_voom"

minCntFlag=c(10)
compFlag3=""

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
        if (colIdPV=="pv") {
            i=which(stat_1$PValue<pThres)
        } else {
            #i=which(stat_1$FDR<pThres)
            i=which(stat_1$Qvalue<pThres)
        }
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
library(marray)
source(paste(dirSrc,"functions/heatmap.5.2.R",sep=""))
source(paste(dirSrc,"functions/heatmapAcgh.7.1.R",sep=""))

outFormat="pdf"
outFormat="png"

sampleBar=""
sampleBar="cluster"

geneBar=""
geneBar="clusterPr"

centrFlag="_noCentering"
centrFlag=""

subsetFlag=""

numPr=500
pThres=10^-8
pThres=10^-6
pThres=0.05

compList="TopVar500"
compList=c("TopVar500","TGFbetaVuntreated_4hrs","TGFbetaVuntreated_8hrs","TGFbetaVuntreated_12hrs","TGFbetaVuntreated_36hrs")

colGeneId="geneId"; colIdPV="FDR"; colNamePV="QV"

tblCC=NULL
for (compFlag in compList) {
    if (length(grep("Rnd",compFlag))==1) {
        rndVec=paste("_rnd",1:4,sep="")
        #		rndVec=paste("_rnd",1:20,sep="")
    } else {
        rndVec=""
    }
    
    for (rndId in rndVec) {
        limFCmmu=c(-6,6)
        if (compFlag%in%c("TGFbetaVuntreated_4hrs","TGFbetaVuntreated_8hrs","TGFbetaVuntreated_12hrs","TGFbetaVuntreated_36hrs")) {
            load(paste(datadir,"dge_Homo_sapiens_",compFlag,".RData",sep=""))
            switch(compFlag,
            "TGFbetaVuntreated_4hrs"={
                stat_1=stat1_4
            },
            "TGFbetaVuntreated_8hrs"={
                stat_1=stat1_8
            },
            "TGFbetaVuntreated_12hrs"={
                stat_1=stat1_12
            },
            "TGFbetaVuntreated_36hrs"={
                stat_1=stat1_36
            }
            )
            timeThis=as.integer(sub("hrs","",strsplit(compFlag,"_")[[1]][2]))
            compName=paste(timeThis,"hr: TGFbeta vs. untreated",sep="")
            cnt_1=t(t(dgeT$counts)*dgeT$samples$norm.factors); phen_1=cbind(phen10[match(rownames(dgeT$samples),phen10$id),],dgeT$samples[,c("lib.size","norm.factors")]); ann_1=ann10[match(rownames(dgeT$counts),ann10$geneId),]
        } else {
            if (substr(compFlag,1,nchar("TopVar"))=="TopVar") {
                compName=paste("Top ",sub("TopVar","",compFlag)," variable genes",sep="")
            }
            phen_1=phen10
            phen_1=phen_1[!phen_1$id%in%c("Donor_1_TGFbeta_4hrs","Donor_3_4_hrs"),]
            stat_1=NULL
            cnt_1=t(t(cnt10[,match(phen_1$id,colnames(cnt10))])*phen_1$norm.factors); ann_1=ann10[match(rownames(cnt10),ann10$geneId),]
        }
        compFlag2=compFlag
        phen_1$id2=phen_1$id
        phen_1$id2=sub("Donor_","",sub("TGFbeta","tgfb",sub("hrs","h",sub("_hrs","hrs",phen_1$id))))
        for (transFlag in c("")) {
            if (transFlag=="") {
                subsetFlag=subsetName=""
            } else {
                subsetFlag=paste("_",tolower(transFlag),sep="")
                subsetName=paste(", ",transFlag,sep="")
            }
            
            if (compFlag%in%c("TGFbetaVuntreated_4hrs","TGFbetaVuntreated_8hrs","TGFbetaVuntreated_12hrs","TGFbetaVuntreated_36hrs")) {
                i1=which(stat_1[,colIdPV]<pThres)
                if (length(i1)==0) next
                fNameOut=paste(compFlag,subsetFlag,centrFlag,"_",colNamePV,pThres,sep="")
                header=paste(compName,subsetName,", ",colNamePV,"<",pThres,sep="")
                dat0=cnt_1[match(stat_1[,colGeneId][i1],rownames(cnt_1)),]
            } else {
                fNameOut=paste(compFlag,subsetFlag,centrFlag,rndId,sep="")
                header=paste(compName,subsetName,sep="")
                dat0=cnt_1
            }
            dat0=dat0+min(c(dat0[dat0!=0]))*0.1
            expr=log2(dat0)
            ann=ann_1[match(rownames(expr),ann_1$geneId),]
            phen=phen_1
            
            i2=1:nrow(expr)
            if (length(grep("Rnd",compFlag))==1) {
                if (length(grep("Rnd500",compFlag))==1) {
                    i1=1:numPr
                }
                header=paste(header,": ",length(i1)," random genes",sep="")
                geneBar="clusterPr"
                #		set.seed(5453)
                i2=sample(1:nrow(expr),length(i1),replace=F)
                expr=expr[i2,]
                ann=ann[i2,]
                i=1:nrow(expr)
            } else if (length(grep("TopVar",compFlag))==1) {
                #header=paste(header,": ",nrow(expr)," top var genes",sep="")
                geneBar="clusterPr"
                varGene=apply(expr,1,var,na.rm=T)
                i2=order(varGene,decreasing=T)[1:numPr]
                expr=expr[i2,]
                ann=ann[i2,]
                i=1:nrow(expr)
            } else if (length(grep("Top",compFlag))==1) {
                header=paste(header,", n=",nrow(expr),sep="")
                geneBar=""
                geneBar="clusterPr"
                i=order(ann$logFC)
            } else {
                geneBar=""
                expr=expr[i2,]
                ann=cbind(ann[i2,],logFC=stat_1$logFC[match(ann[i2,colGeneId],stat_1[,colGeneId])])
                i=order(ann$logFC)
            }
            
            if (transFlag=="") {
                j=1:ncol(expr)
            } else {
                j=which(phen$translocation==transFlag)
            }
            
            
            arrayData=expr[i,j]
            ann=ann[i,]
            phen2=phen[j,]
            
            phenAll=phen10
            phenAll$id2=phenAll$id
            phenAll$id2=sub("Donor_","",sub("TGFbeta","tgfb",sub("hrs","h",sub("_hrs","hrs",phenAll$id))))
            
            if (centrFlag=="") {
                centr=apply(arrayData,1,median,na.rm=T)
                arrayData=arrayData-centr
            }
            
            varList=c("donor","treat","time")
            varName=paste(varList," ",sep="")
            k=which(varList%in%names(phen2))
            varListAll=varList
            varNameAll=varName
            varList=varList[k]
            varName=varName[k]
            
            colList=c("skyblue","blue","yellow","purple","red")
            colList=c("brown","red","orange","yellow","green","cyan","skyblue","blue","pink","magenta","purple","darkgreen")
            colList2=c("skyblue","blue")
            colHM=c("red","blue","grey")
            
            distMethod="pearson"
            linkMethod="ward.D2"
            
            cloneName=ann$geneSymbol
            cloneName=rep("",nrow(ann))
            if (length(grep("Rnd|TopVar",compFlag))==1) {
                cloneCol=NULL
            } else {
                cloneCol=matrix(rep("white",nrow(arrayData)),nrow=1)
                k1=1; kk=which(names(ann)=="logFC")
                x=round(ann[,kk]); x=x-min(x,na.rm=T)+1
                grpUniq=sort(unique(x[!is.na(x)]))
                x=round(ann[,kk]); x[x<limFCmmu[1]]=limFCmmu[1]; x[x>limFCmmu[2]]=limFCmmu[2]; x=x+limFCmmu[2]+1
                grpUniq=limFCmmu[1]:limFCmmu[2]
                cloneColUniq=gray(0:(length(grpUniq)-1)/length(grpUniq))
                cloneCol[k1,]=cloneColUniq[x]
                rownames(cloneCol)="log2FC "
            }
            
            if (subsetFlag=="") {
                samName=rep("",ncol(arrayData))
            } else {
                samName=phen2$id2
            }
            samName=phen2$id2
            samCol=NULL
            samCol=matrix(nrow=length(varList),ncol=nrow(phen2))
            for (varId in 1:length(varList)) {
                if (varList[varId]%in%c("lib.size")) {
                    j=match(phen2$id,phenAll$id)
                    x=round(phenAll[,varList[varId]])
                    lim=range(x,na.rm=T)
                    #lim=quantile(x,probs=c(.1,.9),na.rm=T)
                    x[x<lim[1]]=lim[1]; x[x>lim[2]]=lim[2]
                    grpUniq=lim[1]:lim[2]
                    samColUniq=gray(0:(length(grpUniq)-1)/length(grpUniq))
                    samCol[varId,]=samColUniq[x[j]]
                } else {
                    if (varList[varId]%in%c("time")) {
                        x=phenAll[,varList[varId]]
                    } else {
                        x=as.character(phenAll[,varList[varId]])
                    }
                    x[x==""]=NA; x=as.integer(as.factor(x))
                    grpUniq=sort(unique(x))
                    x=x[match(phen2$id,phenAll$id)]
                    if (length(grpUniq)<=length(colList2)) {
                        samCol[varId,]=colList2[x]
                    } else if (length(grpUniq)<=length(colList)) {
                        samCol[varId,]=colList[x]
                    } else {
                        samCol[varId,]=rainbow(length(grpUniq))[x]
                    }
                }
            }
            rownames(samCol)=varName
            
            print("summary(range(c(arrayData),na.rm=T))")
            print(summary(range(c(arrayData),na.rm=T)))
            if (centrFlag=="") {
                limit=c(-120000,120000)
                limit=c(-10000,10000)
                limit=c(-8,8)
                limit=c(-1,1)
                limit=c(-3,3)
            } else {
                limit=c(8,13)
            }
            main=NULL
            main=header
            
            ncc=ncr=2
            ncc=ncr=NA
            
            switch(distMethod,
            "pearson"={distMat=as.dist(1 - cor(arrayData,method=distMethod,use="complete.obs"))
                if (sampleBar=="cluster") {
                    clustC=hclust(distMat, method=linkMethod)
                } else {
                    clustC=NA
                    ncc=NA
                }
                if (geneBar=="clusterPr") {
                    distMat=as.dist(1 - cor(t(arrayData),method=distMethod,use="complete.obs"))
                    clustR=hclust(distMat, method=linkMethod)
                } else {
                    clustR=NA
                    ncr=NA
                }
            },
            "spearman"={distMat=as.dist(1 - cor(arrayData,method=distMethod,use="complete.obs"))
                distMat=as.dist(1 - cor(t(arrayData),method=distMethod,use="complete.obs"))
            },
            "euclidean"={distMat=dist(t(arrayData), method=distMethod)
                distMat=dist(arrayData, method=distMethod)
            }
            )
            
            if (F) {
                subDir <- paste(compFlag,sep="")
                if (!file.exists(subDir)){
                    dir.create(file.path(subDir))
                }
                subDir=paste(subDir,"/",sep="")
            }
            subDir=""
            if (outFormat=="png") {
                margins=c(6,1)
                margins=c(10,20)
                png(paste(subDir,"heatmap",fNameOut,".png",sep=""),width=480*2,height=480*2)
            } else {
                margins=c(12,5)
                pdf(paste(subDir,"heatmap",fNameOut,".pdf",sep=""))
            }
            #hcc=heatmap3(x=arrayData, Rowv=as.dendrogram(clustR), Colv=as.dendrogram(clustC), distfun=distMethod, hclustfun=hclust, symm=F, ColSideColors=samCol, RowSideColors=cloneCol, labCol=samName, labRow=cloneName, ncr=ncr, ncc=ncc, scale="none", na.rm=F, margins=margins, main=main, xlab=NULL, ylab=NULL, zlm=limit,cexCol=2, , high=colHM[1], low=colHM[2], mid=colHM[3])
            hcc=heatmap3(x=arrayData, Rowv=clustR, Colv=clustC, distfun=distMethod, hclustfun=hclust, symm=F, ColSideColors=samCol, RowSideColors=cloneCol, labCol=samName, labRow=cloneName, ncr=ncr, ncc=ncc, scale="none", na.rm=F, margins=margins, main=main, xlab=NULL, ylab=NULL, zlm=limit,cexCol=2, , high=colHM[1], low=colHM[2], mid=colHM[3])
            dev.off()
        }
    }
}
if (!is.null(samCol)) {
    for (varId in 1:length(varListAll)) {
        if (outFormat=="png") {
            png(paste("heatmapSampleColorBarLegend_",varListAll[varId],".png",sep=""))
        } else {
            pdf(paste("heatmapSampleColorBarLegend_",varListAll[varId],".pdf",sep=""))
        }
        if (varListAll[varId]%in%c("age","wbc")) {
            x=round(phenAll[,varListAll[varId]])
            lim=range(x,na.rm=T)
            #lim=quantile(x,probs=c(.1,.9),na.rm=T)
            grpUniq=lim[1]:lim[2]
            samColUniq=gray(0:(length(grpUniq)-1)/length(grpUniq))
            heatmapColorBar(limit=lim,cols=c(samColUniq[c(length(samColUniq),1)],median(samColUniq)))
        } else {
            if (varList[varId]%in%c("time")) {
                x=phenAll[,varListAll[varId]]
            } else {
                x=as.character(phenAll[,varListAll[varId]]); x[x==""]=NA
            }
            grpUniq=table(x)
            #		grpUniq=paste(names(grpUniq)," (",grpUniq,")",sep="")
            grpUniq=names(grpUniq)
            k=1:length(grpUniq)
            if (length(grpUniq)<=length(colList2)) {
                sampleColorLegend(tls=grpUniq[k],col=colList2,legendTitle=varNameAll[varId])
            } else if (length(grpUniq)<=length(colList)) {
                sampleColorLegend(tls=grpUniq[k],col=colList,legendTitle=varNameAll[varId])
            } else {
                sampleColorLegend(tls=grpUniq[k],col=rainbow(length(grpUniq)),legendTitle=varNameAll[varId])
            }
        }
        dev.off()
    }
}
if (F) {
heatmapColorBar=function(limit,cols=c("green","red","black"),main=NULL) {
    try <- maPalette(high=cols[1], low=cols[2], mid=cols[3],k=15)
    maColorBar(try, scale=limit,k=3,main=main)
}
}
if (outFormat=="png") {
    png(paste("heatmapColorRange.png",sep=""),width=480,height=140)
} else {
    pdf(paste("heatmapColorRange.pdf",sep=""))
}
heatmapColorBar(limit=limit,cols=colHM,main="Heatmap color range")
dev.off()
png("heatmaplog2FoldChangeColorBarRange.png",width=480,height=140)
grpUniq=limFCmmu[1]:limFCmmu[2]
colColUniq=gray(0:(length(grpUniq)-1)/length(grpUniq))
heatmapColorBar(limit=limFCmmu,cols=c(colColUniq[c(length(colColUniq),1)],median(colColUniq)),main="log2FC")
dev.off()

###########################################################
###########################################################


###########################################################
###########################################################