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

## ---------------------------------
#####################################################################################
#####################################################################################


## -------------------
datadir="results/final/misc/"

phen=read.table(paste("docs/sampleInfo.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
ann=read.table(paste(datadir,"ann_Homo_sapiens.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
ann38=ann

## -------------------
library(biomaRt)

dataset <- "hsapiens_gene_ensembl"
ensembl <- useMart("ensembl", dataset = dataset)
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
attrs <- c("ensembl_gene_id", "hgnc_symbol", "hgnc_id", "chromosome_name", "start_position", "end_position","strand", "gene_biotype", "description")
attrs <- c("ensembl_gene_id", "hgnc_symbol", "hgnc_id", "entrezgene", "chromosome_name", "start_position", "end_position","strand", "gene_biotype", "description")
res=getBM(attributes=attrs, filters = "ensembl_gene_id", values = ann$geneId, mart=ensembl, curl = NULL, checkFilters = TRUE, verbose = FALSE, uniqueRows = TRUE, bmHeader = FALSE)
i=match(ann$geneId,res$ensembl_gene_id); i1=which(!is.na(i)); i2=i[i1]; i12=which(is.na(i))
if (length(i12)!=0) {
    tmp=res[1:length(i12),]
    for (k in 1:ncol(tmp)) {
        tmp[,k]=NA
    }
    tmp$ensembl_gene_id=ann$geneId[i12]
    res=rbind(res,tmp)
}
id=unique(res[,1][duplicated(res[,1])])
for (k in 1:length(id)) {
    i=which(res[,1]==id[k])
    if (sum(!duplicated(res$chromosome_name[i]))!=1 | sum(!duplicated(res$start_position[i]))!=1 | sum(!duplicated(res$end_position[i]))!=1) cat(k,"\n")
}
## No replicate gene IDs with different position
i=match(ann$geneId,res$ensembl_gene_id); i1=which(!is.na(i)); i2=i[i1]
res=res[i2,]
ann=res
rm(res)
names(ann)[match(c("ensembl_gene_id"),names(ann))]="geneId"
ann37=ann

## -------------------
ann=ann37

names(ann)[match(c("hgnc_symbol"),names(ann))]=c("geneSym")
ann$chr=as.integer(ann$chromosome_name)
ann$chr[which(ann$chromosome_name=="X")]=23
ann$chr[which(ann$chromosome_name=="Y")]=24
ann2=ann

#####################################################################################
#####################################################################################
cohortFlag="_newman"

## -------------------
datadir=paste("docs/papers/",sub("_","",cohortFlag),"/",sep="")
fName="Newman et al LM22 list nmeth.3337-S2_SuppTable1_DEGs.txt"
candGene=read.table(paste(datadir,fName,sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,skip=2)
names(candGene)[match(c("X22.leukocyte.subsets...activation.states","B.cells.naive","B.cells.memory","Plasma.cells","T.cells.CD8","T.cells.CD4.naive","T.cells.CD4.memory.resting","T.cells.CD4.memory.activated","T.cells.follicular.helper","T.cells.regulatory..Tregs.","T.cells.gamma.delta","NK.cells.resting","NK.cells.activated","Monocytes","Macrophages.M0","Macrophages.M1","Macrophages.M2","Dendritic.cells.resting","Dendritic.cells.activated","Mast.cells.resting","Mast.cells.activated","Eosinophils","Neutrophils"),names(candGene))]=
    c("geneSym","bCellNaive","bCellMem","plasmaCell","tCellCD8","tCellCD4Naive","tCellCD4MemRest","tCellCD4MemActiv","tCellFollHelp","tReg","tCellGammaDelta","nkCellRest","nkCellActiv","monocyte","m0","m1","m2","dendCellRest","dendCellActiv","mastCellRest","mastCellActiv","eosinophil","neutrophil")
table(tbl$monocyte)

#####################################################################################
#####################################################################################
## PLOS paper
## Abbas 2005

## -------------------
if (F) {
    GRCh=34
    build=hg16
    ## Not working
    library(biomaRt)
    ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=GRCh)
}

## -------------------
cohortFlag="_pollara2017"

pThres=0.05

## -------------------
datadir=paste("docs/papers/",sub("_","",cohortFlag),"/",sep="")
fName="journal.pone.0169271.s003_Monocytes.txt"
tbl=read.table(paste(datadir,fName,sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,skip=3)
names(tbl)[match(c("Module.code","Module.description","Module.source","Module.content"),names(tbl))]=c("code","desc","src","content")
for (k in 1:ncol(tbl)) tbl[,k]=gsub("\"","",tbl[,k])
candGene=NULL
colId=which(names(tbl)!="content")
for (k in 1:nrow(tbl)) {
    x=strsplit(gsub(" ","",tbl$content[k]),",")[[1]]
    tbl2=cbind(as.data.frame(matrix(rep(unlist(tbl[k,colId]),length(x)),nrow=length(x),byrow=T),stringsAsFactors=F),geneSym=x)
    candGene=rbind(candGene,tbl2)
}
names(candGene)=c(names(tbl)[colId],"geneSym")

datadir=""
compList=c("antiPd1RespBi")
alt="two.sided"

datadir="results/comparison/"
compList=c("TGFbetaVuntreated_4hrs","TGFbetaVuntreated_8hrs","TGFbetaVuntreated_12hrs","TGFbetaVuntreated_36hrs")
compList=c("TGFbetaVuntreated_4hrs","TGFbetaVuntreated_8hrs","TGFbetaVuntreated_12hrs","TGFbetaVuntreated_36hrs","TGFbetaVuntreated_8hrs_12hrs")
alt="greater"

genesetUniq=unique(candGene$code)
tmp=rep(NA,length(compList)*length(genesetUniq))
tmpC=rep("",length(compList)*length(genesetUniq))
#tbl=data.frame(comp=tmpC,geneset=tmpC,or=tmp,,pv=tmp,stringsAsFactors=F)
orMat=pvMat=matrix(nrow=length(genesetUniq),ncol=length(compList),dimnames=list(genesetUniq,compList))
kk=1
for (cId in 1:length(compList)) {
    compFlag=compList[cId]
    for (gId in 1:length(genesetUniq)) {
        cat("\n\n=========",compFlag,", ",genesetUniq[gId],"=============\n")
        #load(paste(datadir,"dge_Homo_sapiens_",compFlag,".RData",sep=""))
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
        },
        "TGFbetaVuntreated_8hrs_12hrs"={
            stat_1=rbind(stat1_8,stat1_12[!stat1_12$geneId%in%stat1_8$geneId,])
            for (k in c("logFC","logCPM","t")) stat_1[,k]=NA
            for (k in c("PValue","Qvalue")) stat_1[,k]=1
            i11=which(stat1_8$Qvalue<pThres)
            i21=which(stat1_12$Qvalue<pThres)
            i=match(stat1_8$geneId[i11],stat1_12$geneId[i21]); i1=i11[which(!is.na(i))]; i2=i21[i[which(!is.na(i))]]
            i11=i11[!i11%in%i1]
            i21=i21[!i21%in%i2]
            i1=i1[which(sign(stat1_8$logFC[i1])==sign(stat1_12$logFC[i2]))]
            stat_1$Qvalue[which(stat_1$geneId%in%c(stat1_8$geneId[c(i11,i1)],stat1_12$geneId[i21]))]=0
            i1=which(stat_1$Qvalue<pThres & stat_1$geneId%in%stat1_8$geneId[stat1_8$Qvalue<pThres])
            i2=which(stat_1$Qvalue<pThres & stat_1$geneId%in%stat1_12$geneId[stat1_12$Qvalue<pThres])
            i0=which(stat_1$Qvalue>=pThres & stat_1$geneId%in%stat1_8$geneId[stat1_8$Qvalue<pThres] & stat_1$geneId%in%stat1_12$geneId[stat1_12$Qvalue<pThres])
            cat("No. signif 8hr ",length(i1),", 12hr ",length(i2),", signif 8hr/12hr but in opp dir",length(i0),"\n")
        },
        "antiPd1RespBi"={
            stat_1=stat20_1
        }
        )
        nm=data.frame(name1=c("hgnc_symbol","PValue","Qvalue"),name2=c("geneSym","pv","qv"),stringsAsFactors=F); nm=nm[nm$name1%in%names(stat_1),]
        if (nrow(nm)!=0) names(stat_1)[match(nm$name1,names(stat_1))]=nm$name2
        stat_1=stat_1[which(stat_1$geneSym!=""),]
        i=match(tolower(stat_1$geneSym),tolower(candGene$geneSym[which(candGene$code==genesetUniq[gId])])); i1=which(!is.na(i)); i2=which(is.na(i))
        stat_1$candGene=0; stat_1$candGene[i1]=1
        if (F) {
            stat_1$pv[stat_1$candGene==0]=NA
            stat_1$qv=NA
            i=which(!is.na(stat_1$pv))
            if (length(i)!=0) {
                stat_1$qv=qvalue(stat_1$pv)$qvalues
            }
        }
        x=table(signif=stat_1$qv<pThres,candGene=stat_1$candGene)
        print(x)
        if (nrow(x)==2) {
            res=fisher.test(x,alternative=alt)
            pv=res$p.value
            or=res$estimate
        } else {
            pv=or=NA
        }
        if (F) {
            tbl$comp[kk]=compFlag
            tbl$geneset[kk]=genesetUniq[gId]
            tbl$or[kk]=or
            tbl$pv[kk]=pv
            kk=kk+1
        }
        orMat[gId,cId]=or
        pvMat[gId,cId]=pv
        #cat("\nPV ",signif(pv,2),"\n")
        #print(x)
    }
}
## Do Holm's adjustment since the tests are not independent
#library(qvalue)
qvMat=matrix(nrow=length(genesetUniq),ncol=length(compList),dimnames=list(genesetUniq,compList))
for (k in 1:ncol(pvMat)) {
    i=which(!is.na(pvMat[,k]))
    if (length(i)!=0) {
        #qvMat[i,k]=qvalue(pvMat[i,k])$qvalues
        qvMat[i,k]=p.adjust(pvMat[i,k],method="holm")
        #qvMat[i,k]=p.adjust(pvMat[i,k],method="bonf")
    }
}
table(c(qvMat)<pThres)
#round(pvMat,2)
#round(qvMat,2)
