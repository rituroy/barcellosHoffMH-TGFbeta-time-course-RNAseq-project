## Validate our data & GSE78220 data with candidate gene lists

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

if (F) {
    ## -------------------

    phen=read.table(paste("docs/sampleInfo.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)

    source("/Users/royr/Downloads/barcellosHoffMH-TGFbetaTC-tgfb/funcs.R")
    ann38=getAnnotation(build=38)
    ann37=getAnnotation(build=37)

    #####################################################################################
    #####################################################################################

    ## -------------------
    ann=ann37
    ann=ann38
}

#####################################################################################
#####################################################################################
cohortFlag="_newman2015"

## -------------------
datadir=paste("docs/papers/",sub("_","",cohortFlag),"/",sep="")
fName="Newman et al LM22 list nmeth.3337-S2_SuppTable1_DEGs.txt"
candGene=read.table(paste(datadir,fName,sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,skip=2)
names(candGene)[match(c("X22.leukocyte.subsets...activation.states","B.cells.naive","B.cells.memory","Plasma.cells","T.cells.CD8","T.cells.CD4.naive","T.cells.CD4.memory.resting","T.cells.CD4.memory.activated","T.cells.follicular.helper","T.cells.regulatory..Tregs.","T.cells.gamma.delta","NK.cells.resting","NK.cells.activated","Monocytes","Macrophages.M0","Macrophages.M1","Macrophages.M2","Dendritic.cells.resting","Dendritic.cells.activated","Mast.cells.resting","Mast.cells.activated","Eosinophils","Neutrophils"),names(candGene))]=
    c("geneSym","bCellNaive","bCellMem","plasmaCell","tCellCD8","tCellCD4Naive","tCellCD4MemRest","tCellCD4MemActiv","tCellFollHelp","tReg","tCellGammaDelta","nkCellRest","nkCellActiv","monocyte","m0","m1","m2","dendCellRest","dendCellActiv","mastCellRest","mastCellActiv","eosinophil","neutrophil")
table(candGene$monocyte)
x=candGene
colId=2:ncol(x)
candGene=data.frame(geneSym=rep(x$geneSym,length(colId)),group=rep(names(x)[colId],each=nrow(x)),keep=c(as.matrix(x[,colId])),stringsAsFactors=F)
candGene=candGene[which(candGene$keep==1),which(names(candGene)%in%c("geneSym","group"))]
candGeneN=candGene

#####################################################################################
#####################################################################################
cohortFlag="_gentles2015"

## -------------------
datadir=paste("docs/papers/",sub("_","",cohortFlag),"/",sep="")
fName="nm.3909-S4_clusterAssignments.txt"
candGene=read.table(paste(datadir,fName,sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,skip=8)
names(candGene)[match(c("Cluster.ID","Gene.symbol","Compound.score"),names(candGene))]=
c("clustId","geneSym","score")
candGene$group=as.character(candGene$clustId)
candGene$group[which(candGene$clustId==13)]="TGFb"
candGeneG=candGene

#####################################################################################
#####################################################################################
## PLOS paper
## Abbas 2005, GRCh=34, build=hg16


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
names(tbl)[match(c("Module.code","Module.description","Module.source","Module.content"),names(tbl))]=c("group","desc","src","content")
for (k in 1:ncol(tbl)) tbl[,k]=gsub("\"","",tbl[,k])
candGene=NULL
colId=which(names(tbl)!="content")
for (k in 1:nrow(tbl)) {
    x=strsplit(gsub(" ","",tbl$content[k]),",")[[1]]
    tbl2=cbind(as.data.frame(matrix(rep(unlist(tbl[k,colId]),length(x)),nrow=length(x),byrow=T),stringsAsFactors=F),geneSym=x)
    candGene=rbind(candGene,tbl2)
}
names(candGene)=c(names(tbl)[colId],"geneSym")
candGeneP=candGene

#####################################################################################
#####################################################################################
## Correlation with tgf-beta
## Association with response

cohort1Flag="_gonzalezJunca"
cohort1Flag="_GSE78220"
cohort1List=c("_gonzalezJunca","_GSE78220")

cohortFlag="_pollara2017"
cohortFlag="_gentles2015"
cohortFlag="_newman2015"
cohortList=c("_pollara2017","_gentles2015","_newman2015","TGFbetaVuntreated_8hrs16fold_qv0.05")
cohortList=c("_pollara2017","_gentles2015","_newman2015")

cohort1List="_gonzalezJunca"
cohortList="_newman2015"

for (cohort1Flag in cohort1List) {
    switch(cohort1Flag,
        "_gonzalezJunca"={
            datadir="results/comparison/"
            compList=c("TGFbetaVuntreated_4hrs_qv0.05","TGFbetaVuntreated_8hrs_qv0.05","TGFbetaVuntreated_12hrs_qv0.05","TGFbetaVuntreated_36hrs_qv0.05","TGFbetaVuntreated_8hrs_12hrs_qv0.05_sameDir")
            compList=c("TGFbetaVuntreated_8hrs16fold_qv0.05","TGFbetaVuntreated_12hrs8fold_qv0.05","TGFbetaVuntreated_8hrs16fold_12hrs8fold_qv0.05_sameDir")
            compList=c("TGFbetaVuntreated_8hrs16fold_qv0.05","TGFbetaVuntreated_12hrs8fold_qv0.05","TGFbetaVuntreated_8hrs16fold_12hrs8fold_qv0.05_sameDir","TGFbetaVuntreated_8hrs_qv0.05","TGFbetaVuntreated_12hrs_qv0.05","TGFbetaVuntreated_8hrs_12hrs_qv0.05_sameDir")
            compList=c("TGFbetaVuntreated_8hrs_qv0.05")
            alt="greater"
            alt="two.sided"
            multtestFlag="_allFeatures"
        },
        "_GSE78220"={
            datadir=""
            compList=c("respVprogDisease_pv0.05")
            alt="two.sided"
            multtestFlag="_allFeatures"
            multtestFlag="_candGene"
            multtestFlag="_none"
        }
    )
    for (cohortFlag in cohortList) {
        if (cohort1Flag=="_gonzalezJunca" & substr(cohortFlag,1,nchar("TGFbetaVuntreated_"))=="TGFbetaVuntreated_") next
        cat("\n\n=========",cohort1Flag,", ",cohortFlag,"=============\n",sep="")
        switch(cohortFlag,
            "_pollara2017"={
                candGene=candGeneP
                genesetUniq=unique(candGene$group)
            },
            "_gentles2015"={
                candGene=candGeneG
                genesetUniq="TGFb"
            },
            "_newman2015"={
                candGene=candGeneN
                genesetUniq=unique(candGene$group)
            },
            "TGFbetaVuntreated_8hrs16fold_qv0.05"={
                stat_1=stat1_8_10
                for (k in c("PValue","Qvalue")) stat_1[,k]=1
                i=which(stat1_8_10$Qvalue<pThres & abs(stat1_8_10$logFC)>=4)
                stat_1$Qvalue[i]=0
                nm=data.frame(name1=c("hgnc_symbol","PValue","Qvalue"),name2=c("geneSym","pv","qv"),stringsAsFactors=F); nm=nm[nm$name1%in%names(stat_1),]
                if (nrow(nm)!=0) names(stat_1)[match(nm$name1,names(stat_1))]=nm$name2
                candGene=stat_1
                candGene$group=""
                candGene$group[i]="TGFb"
                genesetUniq="TGFb"
            }
        )

        tmp=rep(NA,length(compList)*length(genesetUniq))
        tmpC=rep("",length(compList)*length(genesetUniq))
        #tbl=data.frame(comp=tmpC,geneset=tmpC,or=tmp,lcbtmp,ucb=tmp,pv=tmp,stringsAsFactors=F)
        orMat=pvMat=matrix(nrow=length(genesetUniq),ncol=length(compList),dimnames=list(genesetUniq,compList))
        kk=1
        tbl1=NULL
        for (cId in 1:length(compList)) {
            compFlag=compList[cId]
            for (gId in 1:length(genesetUniq)) {
                #cat("\n\n=========",cohort1Flag,", ",compFlag,", ",cohortFlag,", ",genesetUniq[gId],"=============\n",sep="")
                #load(paste(datadir,"dge_Homo_sapiens_",compFlag,".RData",sep=""))
                switch(compFlag,
                "TGFbetaVuntreated_4hrs_qv0.05"={
                    stat_1=stat1_4_10
                },
                "TGFbetaVuntreated_8hrs_qv0.05"={
                    stat_1=stat1_8_10
                },
                "TGFbetaVuntreated_12hrs_qv0.05"={
                    stat_1=stat1_12_10
                },
                "TGFbetaVuntreated_36hrs_qv0.05"={
                    stat_1=stat1_36_10
                },
                "TGFbetaVuntreated_8hrs_12hrs_qv0.05_sameDir"={
                    stat_1=rbind(stat1_8_10,stat1_12_10[!stat1_12_10$geneId%in%stat1_8_10$geneId,])
                    for (k in c("logFC","logCPM","t")) stat_1[,k]=NA
                    for (k in c("PValue","Qvalue")) stat_1[,k]=1
                    i11=which(stat1_8_10$Qvalue<pThres)
                    i21=which(stat1_12_10$Qvalue<pThres)
                    i=match(stat1_8_10$geneId[i11],stat1_12_10$geneId[i21]); i1=i11[which(!is.na(i))]; i2=i21[i[which(!is.na(i))]]
                    i11=i11[!i11%in%i1]
                    i21=i21[!i21%in%i2]
                    i1=i1[which(sign(stat1_8_10$logFC[i1])==sign(stat1_12_10$logFC[i2]))]
                    stat_1$Qvalue[which(stat_1$geneId%in%c(stat1_8_10$geneId[c(i11,i1)],stat1_12_10$geneId[i21]))]=0
                    i1=which(stat_1$Qvalue<pThres & stat_1$geneId%in%stat1_8_10$geneId[stat1_8_10$Qvalue<pThres])
                    i2=which(stat_1$Qvalue<pThres & stat_1$geneId%in%stat1_12_10$geneId[stat1_12_10$Qvalue<pThres])
                    i0=which(stat_1$Qvalue>=pThres & stat_1$geneId%in%stat1_8_10$geneId[stat1_8_10$Qvalue<pThres] & stat_1$geneId%in%stat1_12_10$geneId[stat1_12_10$Qvalue<pThres])
                    #cat("No. signif 8hr ",length(i1),", 12hr ",length(i2),", signif 8hr/12hr but in opp dir",length(i0),"\n")
                },
                "TGFbetaVuntreated_8hrs16fold_qv0.05"={
                    stat_1=stat1_8_10
                    for (k in c("PValue","Qvalue")) stat_1[,k]=1
                    i=which(stat1_8_10$Qvalue<pThres & abs(stat1_8_10$logFC)>=4)
                    stat_1$Qvalue[i]=0
                },
                "TGFbetaVuntreated_12hrs8fold_qv0.05"={
                    stat_1=stat1_12_10
                    for (k in c("PValue","Qvalue")) stat_1[,k]=1
                    i=which(stat1_12_10$Qvalue<pThres & abs(stat1_12_10$logFC)>=3)
                    stat_1$Qvalue[i]=0
                },
                "TGFbetaVuntreated_8hrs16fold_12hrs8fold_qv0.05_sameDir"={
                    stat_1=rbind(stat1_8_10,stat1_12_10[!stat1_12_10$geneId%in%stat1_8_10$geneId,])
                    for (k in c("logFC","logCPM","t")) stat_1[,k]=NA
                    for (k in c("PValue","Qvalue")) stat_1[,k]=1
                    i11=which(stat1_8_10$Qvalue<pThres & abs(stat1_8_10$logFC)>=4)
                    i21=which(stat1_12_10$Qvalue<pThres & abs(stat1_12_10$logFC)>=3)
                    i=match(stat1_8_10$geneId[i11],stat1_12_10$geneId[i21]); i1=i11[which(!is.na(i))]; i2=i21[i[which(!is.na(i))]]
                    i11=i11[!i11%in%i1]
                    i21=i21[!i21%in%i2]
                    i1=i1[which(sign(stat1_8_10$logFC[i1])==sign(stat1_12_10$logFC[i2]))]
                    stat_1$Qvalue[which(stat_1$geneId%in%c(stat1_8_10$geneId[c(i11,i1)],stat1_12_10$geneId[i21]))]=0
                    i1=which(stat_1$Qvalue<pThres & stat_1$geneId%in%stat1_8_10$geneId[stat1_8_10$Qvalue<pThres & abs(stat1_8_10$logFC)>=4])
                    i2=which(stat_1$Qvalue<pThres & stat_1$geneId%in%stat1_12_10$geneId[stat1_12_10$Qvalue<pThres & abs(stat1_12_10$logFC)>=3])
                    i0=which(stat_1$Qvalue>=pThres & stat_1$geneId%in%stat1_8_10$geneId[stat1_8_10$Qvalue<pThres & abs(stat1_8_10$logFC)>=4] & stat_1$geneId%in%stat1_12_10$geneId[stat1_12_10$Qvalue<pThres & abs(stat1_12_10$logFC)>=3])
                    #cat("No. signif 8hr ",length(i1),", 12hr ",length(i2),", signif 8hr/12hr but in opp dir",length(i0),"\n")
                },
                "respVprogDisease_pv0.05"={
                    stat_1=stat20_1
                }
                )
                nm=data.frame(name1=c("hgnc_symbol","PValue","Qvalue"),name2=c("geneSym","pv","qv"),stringsAsFactors=F); nm=nm[nm$name1%in%names(stat_1),]
                if (nrow(nm)!=0) names(stat_1)[match(nm$name1,names(stat_1))]=nm$name2
                stat_1=stat_1[which(stat_1$geneSym!="" & !is.na(stat_1$pv)),]
                i=match(tolower(stat_1$geneSym),tolower(candGene$geneSym[which(candGene$group==genesetUniq[gId])])); i1=which(!is.na(i)); i2=which(is.na(i))
                stat_1$candGene=0; stat_1$candGene[i1]=1
                if (multtestFlag=="_candGene") {
                    stat_1$pv[stat_1$candGene==0]=NA
                    stat_1$qv=NA
                    i=which(!is.na(stat_1$pv))
                    if (length(i)!=0) {
                        stat_1$qv=qvalue(stat_1$pv)$qvalues
                    }
                } else if (multtestFlag=="_none") {
                    stat_1$qv=stat_1$pv
                }
                x=table(signif=stat_1$qv<pThres,candGene=stat_1$candGene)
                #print(x)
                numFeat=nrow(stat_1)
                numSignif=sum(stat_1$qv<pThres)
                numInSignat=sum(stat_1$candGene==1)
                genesInSignat=paste(stat_1$geneSym[which(stat_1$candGene==1)],collapse="/")
                numSignifAndInSignat=sum(stat_1$qv<pThres & stat_1$candGene==1)
                if (nrow(x)==2) {
                    res=fisher.test(x,alternative=alt)
                    pv=res$p.value
                    or=res$estimate
                    ci=res$conf.int
                } else {
                    pv=or=ci=NA
                }
                if (F) {
                    tbl$comp[kk]=compFlag
                    tbl$geneset[kk]=genesetUniq[gId]
                    tbl$or[kk]=or
                    tbl$lcb[kk]=ci[1]
                    tbl$ucb[kk]=ci[2]
                    tbl$pv[kk]=pv
                    kk=kk+1
                }
                orMat[gId,cId]=or
                pvMat[gId,cId]=pv
                #cat("\nPV ",signif(pv,2),"\n")
                #print(x)
                tbl2=c(sub("_","",cohortFlag),genesetUniq[gId],sub("_","",cohort1Flag),compFlag,sub("_","",multtestFlag),numFeat,numSignif,numInSignat,genesInSignat,numSignifAndInSignat,log(or),log(ci),pv,alt)
                tbl1=rbind(tbl1,tbl2)
            }
            grpInfo=data.frame(grp=c("monocyte","m0","neutrophil","dendCellRest","dendCellActiv","m2"),name=c("Monocytes","Macrophages M0","Neutrophils","Dendritic cells resting","Dendritic cells activated","MacrophagesM2"),name2=c("Monocytes","Macrophages\nM0","Neutrophils","Dendritic cells\nresting","Dendritic cells\nactivated","Macrophages\nM2"),stringsAsFactors=F)
            k=which(paste("_",tbl1[,1],sep="")==cohortFlag & paste("_",tbl1[,3],sep="")==cohort1Flag)
            k1=match(grpInfo$grp,tbl1[k,2]); k1=k1[which(!is.na(k1))]
            #k1=1:length(k)
            if (length(k1)!=0) {
                k=k[k1]
                x=apply(tbl1[k,11:13],c(1,2),as.numeric)
                xlim=c(0.5,length(k)+0.5)
                ylim=range(x[is.finite(x)],na.rm=T)
                png(paste("lorPlot_",compFlag,"_val",cohort1Flag,"_signat",cohortFlag,".png",sep=""))
                par(mar=c(5, 4, 4, 2) + 0.1)
                par(mar=c(8, 5, 1, 1) + 0.1)
                plot(1:length(k),x[,1],xlim=xlim,ylim=ylim,cex=2,xaxt="n",xlab="",ylab="log (odds ratio)")
                axis(side=1,at=1:length(k),labels=grpInfo$name2[match(tbl1[k,2],grpInfo$grp)],las=3)
                ci=x[,2:3]; ci[which(ci==-Inf)]=lim[1]-1; ci[which(ci==Inf)]=lim[2]+1
                for (k1 in 1:length(k)) {
                    lines(rep(k1,2),y=ci[k1,],lty="solid")
                    lines(k1+c(-0.1,0.1),y=rep(x[k1,2],2),lty="solid")
                    lines(k1+c(-0.1,0.1),y=rep(x[k1,3],2),lty="solid")
                }
                abline(h=0,lty="dotted")
                dev.off()
            }
        }
        rownames(tbl1)=NULL
        tbl1=as.data.frame(tbl1,stringsAsFactors=F)
        names(tbl1)=c("cohortSignat","geneset","cohortVal","assoWith","multTest","numFeat","numSignif","numSignat","genesSignat","numSignifSignat","lor","lcb","ucb","pv","testAlt")
        for (k in c("numFeat","numSignif","numSignat","numSignifSignat","lor","lcb","ucb","pv")) tbl1[,k]=as.numeric(tbl1[,k])

        ## Do Holm's adjustment since the tests are not independent
        #library(qvalue)
        tbl1$qv=NA
        qvMat=matrix(nrow=length(genesetUniq),ncol=length(compList),dimnames=list(genesetUniq,compList))
        for (k in 1:ncol(pvMat)) {
            i=which(!is.na(pvMat[,k]))
            if (length(i)!=0) {
                #qvMat[i,k]=qvalue(pvMat[i,k])$qvalues
                qvMat[i,k]=p.adjust(pvMat[i,k],method="holm")
                #qvMat[i,k]=p.adjust(pvMat[i,k],method="bonf")
                tbl1$qv[which(tbl1$geneset==rownames(qvMat)[i] & tbl1$assoWith==colnames(qvMat)[k])]=qvMat[i,k]
            }
        }
        if (T) {
            #cat("\ntable(c(qvMat)<pThres)")
            #print(table(c(qvMat)<pThres))
            #round(pvMat,2)
            #round(qvMat,2)
            for (k in 1:ncol(qvMat)) {
                x=table(qvMat[,k]<pThres)
                names(x)=paste(ifelse(multtestFlag=="_none","PV","QV"),c(">=","<"),pThres,sep="")[match(names(x),c("FALSE","TRUE"))]
                cat("\n",colnames(qvMat)[k],":\n",sep="")
                print(x)
            }
        }
        if (F) {
            k=1
            x=as.integer(qvMat[,k]<pThres)
            cat("\nSum ",sum(x),":\n",sep="")
            for (k in 2:ncol(qvMat)) {
                cat(sum(x==1 & qvMat[,k]<pThres),"\n",sep="")
                i=which(x==1 & qvMat[,k]>=pThres)
            }
        }
        names(tbl1)[match(c("multTest","qv"),names(tbl1))]=c("multtestInCohortVal","fdr")
        write.table(tbl1,file=paste("stat_val",cohort1Flag,"_signat",cohortFlag,".txt",sep=""),col.names=T,row.names=F,sep="\t",quote=F)
    }
}

"

=========_gonzalezJunca, _pollara2017=============

TGFbetaVuntreated_8hrs16fold_qv0.05:
QV>=0.05  QV<0.05
9        1

TGFbetaVuntreated_12hrs8fold_qv0.05:
QV>=0.05
10

TGFbetaVuntreated_8hrs16fold_12hrs8fold_qv0.05_sameDir:
QV>=0.05
10

TGFbetaVuntreated_8hrs_qv0.05:
QV>=0.05  QV<0.05
1        9

TGFbetaVuntreated_12hrs_qv0.05:
QV>=0.05  QV<0.05
1        9

TGFbetaVuntreated_8hrs_12hrs_qv0.05_sameDir:
QV>=0.05  QV<0.05
1        9


=========_gonzalezJunca, _gentles2015=============

TGFbetaVuntreated_8hrs16fold_qv0.05:
QV>=0.05
1

TGFbetaVuntreated_12hrs8fold_qv0.05:
QV<0.05
1

TGFbetaVuntreated_8hrs16fold_12hrs8fold_qv0.05_sameDir:
QV<0.05
1

TGFbetaVuntreated_8hrs_qv0.05:
QV<0.05
1

TGFbetaVuntreated_12hrs_qv0.05:
QV<0.05
1

TGFbetaVuntreated_8hrs_12hrs_qv0.05_sameDir:
QV<0.05
1


=========_gonzalezJunca, _newman2015=============

TGFbetaVuntreated_8hrs16fold_qv0.05:
QV>=0.05
22

TGFbetaVuntreated_12hrs8fold_qv0.05:
QV>=0.05
22

TGFbetaVuntreated_8hrs16fold_12hrs8fold_qv0.05_sameDir:
QV>=0.05
22

TGFbetaVuntreated_8hrs_qv0.05:
QV>=0.05  QV<0.05
11       11

TGFbetaVuntreated_12hrs_qv0.05:
QV>=0.05  QV<0.05
11       11

TGFbetaVuntreated_8hrs_12hrs_qv0.05_sameDir:
QV>=0.05  QV<0.05
11       11


===============================================
===============================================
===============================================

=========_GSE78220, _pollara2017=============

respVprogDisease_pv0.05:
PV>=0.05
10


=========_GSE78220, _gentles2015=============

respVprogDisease_pv0.05:
PV<0.05
1


=========_GSE78220, _newman2015=============

respVprogDisease_pv0.05:
PV>=0.05
22

"
