###########################################################
###########################################################

library(marray)
#source(paste(dirSrc,"functions/heatmap.5.2.R",sep=""))
#source(paste(dirSrc,"functions/heatmapAcgh.7.1.R",sep=""))
source(paste(dirSrc,"functions/heatmap.5.6.R",sep=""))
source(paste(dirSrc,"functions/heatmapAcgh.7.3.R",sep=""))

datadirG=""

outFormat="pdf"
outFormat="png"

sampleBar=""
sampleBar="cluster"

geneBar="clusterPr"
geneBar=""

centrFlag="_noCentering"
centrFlag=""
scaleFlag="_scale"
scaleFlag=""

subsetFlag=""

numPr=500
pThres=10^-8
pThres=10^-6
pThres=0.05

limSc=c(-100,100)
limFCmmu=c(-8,8)

cohortFlag=""
compList="TopVar500"
compList=c("TopVar500","TGFbetaVuntreated_4hrs","TGFbetaVuntreated_8hrs","TGFbetaVuntreated_12hrs","TGFbetaVuntreated_36hrs")
subsetList=""

cohortFlag="_GSE33331"
compList=""
#geneList=c("TGFbetaVuntreated_8hrs_qv0.05","TGFbetaVuntreated_8hrs8fold_qv0.05","TGFbetaVuntreated_8hrs16fold_qv0.05","TGFbetaVuntreated_8hrs_12hrs_qv0.05_sameDir","TGFbetaVuntreated_8hrs2fold_12hrs_qv0.05_sameDir","TGFbetaVuntreated_8hrs2fold_12hrs2fold_qv0.05_sameDir","TGFbetaVuntreated_8hrs8fold_12hrs8fold_qv0.05_sameDir")
geneFlag="TGFbetaVuntreated_8hrs16fold_qv0.05"
subsetList=""

cohortFlag="_GSE78220"
compList="RndGene"
compList=""
#geneList=c("TGFbetaVuntreated_8hrs_qv0.05","TGFbetaVuntreated_8hrs8fold_qv0.05","TGFbetaVuntreated_8hrs16fold_qv0.05","TGFbetaVuntreated_8hrs_12hrs_qv0.05_sameDir","TGFbetaVuntreated_8hrs2fold_12hrs_qv0.05_sameDir","TGFbetaVuntreated_8hrs2fold_12hrs2fold_qv0.05_sameDir","TGFbetaVuntreated_8hrs8fold_12hrs8fold_qv0.05_sameDir")
subsetFlag="_ucla"
subsetFlag=""
subsetFlag="_biopsyPreTreat"
subsetFlag="_ucla_biopsyPreTreat"
subsetList=c("","_ucla","_biopsyPreTreat","_ucla_biopsyPreTreat")
geneFlag="TGFbetaVuntreated_8hrs16fold_qv0.05"
subsetList="_ucla_biopsyPreTreat_noOutlierScore"
subsetList="_ucla_biopsyPreTreat"
geneFlag="TGFbetaVuntreated_12hrs8fold_qv0.05"
geneFlag="TGFbetaVuntreated_8hrs16fold_qv0.05"
subsetList="_biopsyPreTreat"

switch(cohortFlag,
"_GSE33331"={
    datV=datObj31$dat
    phenV=datObj31$phen
    annV=datObj31$ann
    scoreMat=datObj31$score
},
"_GSE78220"={
    datV=datObj20$dat
    phenV=datObj20$phen
    annV=datObj20$ann
    scoreMat=datObj20$score
}
)

colGeneId="geneId"; colIdPV="FDR"; colNamePV="QV"

for (geneFlag in c("TGFbetaVuntreated_12hrs8fold_qv0.05","TGFbetaVuntreated_8hrs16fold_qv0.05")) {

tblCC=NULL
for (subsetFlag in subsetList) {
    for (compFlag in compList) {
        if (length(grep("Rnd",compFlag))==1) {
            rndVec=paste("_rnd",1:4,sep="")
            #		rndVec=paste("_rnd",1:20,sep="")
        } else {
            rndVec=""
        }
        
        for (rndId in rndVec) {
            if (cohortFlag%in%c("_GSE78220","_GSE33331")) {
                
            } else {
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
            }
            for (transFlag in c("")) {
                if (cohortFlag%in%c("_GSE78220","_GSE33331")) {
                    subsetName=""
                    tbl=read.table(paste(datadirG,"geneList",cohortFlag,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
                    tbl=tbl[which(tbl$geneList==geneFlag),]
                    nm="weight"; nm1=nm
                    nm="logFC"; nm1="log2FC"
                    tbl=tbl[order(tbl[,nm]),]
                    iC=match(tbl$geneSym,annV$geneSym)
                    if (cohortFlag%in%c("_GSE78220")) {
                        j=which(!duplicated(phenV$patientId))
                        if (subsetFlag!="") {
                            subsetName=paste(", ",subsetFlag,sep="")
                            if (subsetFlag=="_biopsyPreTreat") {
                                j=j[which(phenV$biopsyTime[j]=="pre-treatment")]
                            }
                            if (subsetFlag=="_ucla") {
                                j=j[which(phenV$studySite[j]=="UCLA")]
                            }
                            if (subsetFlag=="_ucla_biopsyPreTreat") {
                                j=j[which(phenV$studySite[j]=="UCLA")]
                                j=j[which(phenV$biopsyTime[j]=="pre-treatment")]
                            }
                            if (subsetFlag=="_ucla_biopsyPreTreat_noOutlierScore") {
                                j=j[which(phenV$studySite[j]=="UCLA")]
                                j=j[which(phenV$biopsyTime[j]=="pre-treatment")]
                                tbl=read.table(paste(datadirG,"scoreMat",cohortFlag,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
                                nm1=paste("t_",geneFlag,sep="")
                                tbl=tbl[match(phenV$id,tbl$id),which(names(tbl)==nm1)]
                                j=j[which(abs(tbl[j])<=90)]
                            }
                        }
                        varList=c("antiPd1Resp")
                        varList=c("antiPd1Resp","studySite","gender","diseaseStatus","vitalStatus","previousMapki","mutation","biopsyTime")
                        x=datV; x[datV==0]=NA
                        x=min(c(x),na.rm=T)/10
                    } else if (cohortFlag%in%c("_GSE33331")) {
                        j=1:nrow(phenV)
                        varList=c("os")
                        x=0
                    } else {
                        j=NULL
                    }
                    if (length(grep("Rnd",compFlag))==1) {
                        #set.seed(5453)
                        iC=sample(1:nrow(datV),length(iC),replace=F)
                        header=paste(length(iC)," random probesets: ",sub("_","",cohortFlag),subsetFlag,sep="")
                        fNameOut=paste("_",geneFlag,compFlag,cohortFlag,subsetFlag,sep="")
                        fNameOut=paste(fNameOut,rndId,sep="")
                        geneBar="clusterPr"
                    } else {
                        header=paste(sub("TGFbetaVuntreated","TGFbVuntreat",geneFlag),": ",sub("_","",cohortFlag),subsetFlag,sep="")
                        fNameOut=paste("_",geneFlag,compFlag,cohortFlag,subsetFlag,sep="")
                    }
                    arrayData=log2(datV[iC,j]+x)
                    annRow=cbind(annV[iC,],tbl[,nm])
                    names(annRow)=c(names(annV),nm1)
                    annRowAll=annRow
                    annCol=phenV[j,]
                    annColAll=phenV
                    #names(annCol)[grep("t_",names(annCol))]="score"
                    #names(annColAll)[grep("t_",names(annColAll))]="score"
                    varName=paste(varList," ",sep="")
                    if (cohortFlag%in%c("_GSE78220")) {
                        tbl=read.table(paste(datadirG,"scoreMat",cohortFlag,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
                        #tbl=tbl[match(phenV$id,tbl$id),which(substr(names(tbl),1,nchar("t_"))=="t_")]
                        nm1=paste("t_",geneFlag,sep="")
                        tbl=tbl[match(phenV$id,tbl$id),which(names(tbl)==nm1)]
                        nm1="score"
                        if (class(tbl)=="numeric") {
                            nm=c(names(annColAll),nm1)
                            annColAll=cbind(annColAll,tbl)
                            names(annColAll)=nm
                            nm=c(names(annCol),nm1)
                            annCol=cbind(annCol,nm=tbl[j])
                            names(annCol)=nm
                        } else {
                            annColAll=cbind(annColAll,tbl)
                            annCol=cbind(annCol,tbl[j,])
                        }
                        varList=c(varList,nm1)
                        varName=c(varName,"score ")
                    }
                } else {
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
                    annRow=ann[i,]
                    annCol=phen[j,]
                    
                    annColAll=phen10
                    annColAll$id2=annColAll$id
                    annColAll$id2=sub("Donor_","",sub("TGFbeta","tgfb",sub("hrs","h",sub("_hrs","hrs",annColAll$id))))
                    
                    varList=c("donor","treat","time")
                    varName=paste(varList," ",sep="")
                }
                
                if (centrFlag=="") {
                    centr=apply(arrayData,1,median,na.rm=T)
                    arrayData=arrayData-centr
                }
                if (scaleFlag=="_scale") {
                    centr=apply(arrayData,1,mad,na.rm=T)
                    arrayData=arrayData/centr
                }
                arrayData[is.infinite(arrayData)]=NA
                i=apply(arrayData,1,function(x) mean(!is.na(x) & !is.infinite(x))!=0)
                arrayData=arrayData[i,]
                annRow=annRow[i,]
                
                k=which(varList%in%names(annCol))
                varListAll=varList
                varNameAll=varName
                varList=varList[k]
                varName=varName[k]
                varFList=varFListAll="weight"
                varFList=varFListAll="log2FC"
                varFName=varFNameAll=paste(varFList," ",sep="")
                
                colList=c("skyblue","blue","yellow","purple","red")
                colList=c("brown","red","orange","yellow","green","cyan","skyblue","blue","pink","magenta","purple","darkgreen")
                colList2=c("skyblue","blue")
                colHM=c("red","blue","grey")
                
                distMethod="pearson"
                linkMethod="ward.D2"
                
                if (cohortFlag%in%c("_GSE78220","_GSE33331")) {
                    distMethod="euclidean"
                    if (cohortFlag%in%c("_GSE78220")) {
                        distMethod="cosine"
                        distMethod="pearson"
                        distMethod="euclidean"
                    }
                    cloneName=rep("",nrow(annRow))
                    cloneName=annRow$geneSym
                    #if (geneBar=="clusterPr") {
                    #    cloneCol=NULL
                    #} else {
                    for (varId in 1:length(varFList)) {
                        cloneCol=matrix(rep("white",nrow(arrayData)),nrow=1)
                        k1=1; kk=which(names(annRow)==varFList[varId])
                        x=round(annRow[,kk])
                        lim=max(abs(x),na.rm=T); lim=c(-lim,lim)
                        lim=lim-min(x,na.rm=T)+1; x=x-min(x,na.rm=T)+1
                        grpUniq=sort(unique(x[!is.na(x)]))
                        lim=limFCmmu
                        x=round(annRow[,kk]); x[x<lim[1]]=lim[1]; x[x>lim[2]]=lim[2]; x=x+lim[2]+1
                        grpUniq=lim[1]:lim[2]
                        cloneColUniq=gray(0:(length(grpUniq)-1)/length(grpUniq))
                        cloneCol[k1,]=cloneColUniq[x]
                    }
                    rownames(cloneCol)=varFName
                        #}
                } else {
                    cloneName=annRow$geneSymbol
                    cloneName=rep("",nrow(annRow))
                    if (length(grep("Rnd|TopVar",compFlag))==1) {
                        cloneCol=NULL
                    } else {
                        cloneCol=matrix(rep("white",nrow(arrayData)),nrow=1)
                        k1=1; kk=which(names(annRow)=="logFC")
                        x=round(annRow[,kk]); x=x-min(x,na.rm=T)+1
                        grpUniq=sort(unique(x[!is.na(x)]))
                        x=round(annRow[,kk]); x[x<limFCmmu[1]]=limFCmmu[1]; x[x>limFCmmu[2]]=limFCmmu[2]; x=x+limFCmmu[2]+1
                        grpUniq=limFCmmu[1]:limFCmmu[2]
                        cloneColUniq=gray(0:(length(grpUniq)-1)/length(grpUniq))
                        cloneCol[k1,]=cloneColUniq[x]
                        rownames(cloneCol)="log2FC "
                    }
                }
                
                if (subsetFlag=="") {
                    samName=rep("",ncol(arrayData))
                } else {
                    samName=annCol$id
                }
                samName=annCol$id
                samCol=NULL
                samCol=matrix(nrow=length(varList),ncol=nrow(annCol))
                for (varId in 1:length(varList)) {
                    if (varList[varId]%in%c("os") | length(grep("t_|score",varList[varId]))==1) {
                        j=match(annCol$id,annColAll$id)
                        x=round(annColAll[,varList[varId]])
                        if (length(grep("t_|score",varList[varId]))==1) {
                            lim=limSc+101
                            x=x+101
                        } else {
                            lim=range(x,na.rm=T)
                            #lim=quantile(x,probs=c(.1,.9),na.rm=T)
                        }
                        x[x<lim[1]]=lim[1]; x[x>lim[2]]=lim[2]
                        grpUniq=lim[1]:lim[2]
                        samColUniq=gray(0:(length(grpUniq)-1)/length(grpUniq))
                        samCol[varId,]=samColUniq[x[j]]
                    } else {
                        if (varList[varId]%in%c("time")) {
                            x=annColAll[,varList[varId]]
                        } else {
                            x=as.character(annColAll[,varList[varId]])
                        }
                        x[x==""]=NA; x=as.integer(as.factor(x))
                        grpUniq=sort(unique(x))
                        x=x[match(annCol$id,annColAll$id)]
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
                if (cohortFlag%in%c("_GSE78220","_GSE33331")) {
                    limit=c(-1,1)
                } else {
                    if (centrFlag=="") {
                        limit=c(-120000,120000)
                        limit=c(-10000,10000)
                        limit=c(-8,8)
                        limit=c(-1,1)
                        limit=c(-3,3)
                    } else {
                        limit=c(8,13)
                    }
                }
                main=NULL
                main=header
                
                nClust=c(NA,NA)
                nClust=c(2,2)
                
                if (sampleBar=="cluster") {
                    switch(distMethod,
                    "pearson"={
                        distMat=as.dist(1 - cor(arrayData,method=distMethod,use="complete.obs"))
                    },
                    "spearman"={
                        distMat=as.dist(1 - cor(arrayData,method=distMethod,use="complete.obs"))
                    },
                    "euclidean"={
                        distMat=dist(t(arrayData), method=distMethod)
                    },
                    "cosine"={
                        if (any(apply(arrayData,2,sum,na.rm=T)==0)) arrayData=arrayData+1
                        distMat=getCosineDist(t(arrayData))
                    }
                    )
                    clustC=hclust(distMat, method=linkMethod)
                } else {
                    clustC=NA
                    nClust[2]=NA
                }
                if (geneBar=="clusterPr") {
                    switch(distMethod,
                    "pearson"={
                        distMat=as.dist(1 - cor(t(arrayData),method=distMethod,use="complete.obs"))
                    },
                    "spearman"={
                        distMat=as.dist(1 - cor(t(arrayData),method=distMethod,use="complete.obs"))
                    },
                    "euclidean"={
                        distMat=dist(arrayData, method=distMethod)
                    },
                    "cosine"={
                        if (any(apply(arrayData,2,sum,na.rm=T)==0)) arrayData=arrayData+1
                        distMat=getCosineDist(arrayData)
                    }
                    )
                    clustR=hclust(distMat, method=linkMethod)
                } else {
                    clustR=NA
                    nClust[1]=NA
                }
                
                if (F) {
                    subDir <- paste(compFlag,sep="")
                    if (!file.exists(subDir)){
                        dir.create(file.path(subDir))
                    }
                    subDir=paste(subDir,"/",sep="")
                }
                subDir=""
                if (outFormat=="png") {
                    if (cohortFlag%in%c("_GSE78220","_GSE33331")) {
                        margins=c(6,1)
                    } else {
                        margins=c(6,1)
                        margins=c(10,20)
                    }
                    png(paste(subDir,"heatmap",fNameOut,".png",sep=""),width=480*2,height=480*2)
                } else {
                    margins=c(12,5)
                    pdf(paste(subDir,"heatmap",fNameOut,".pdf",sep=""))
                }
                #hcc=heatmap3(x=arrayData, Rowv=as.dendrogram(clustR), Colv=as.dendrogram(clustC), distfun=distMethod, hclustfun=hclust, symm=F, ColSideColors=samCol, RowSideColors=cloneCol, labCol=samName, labRow=cloneName, ncr=nClust[1], ncc=nClust[2], scale="none", na.rm=F, margins=margins, main=main, xlab=NULL, ylab=NULL, zlm=limit,cexCol=2, high=colHM[1], low=colHM[2], mid=colHM[3])
                hcc=heatmap3(x=arrayData, Rowv=clustR, Colv=clustC, distfun=distMethod, hclustfun=hclust, symm=F, ColSideColors=samCol, RowSideColors=cloneCol, labCol=samName, labRow=cloneName, ncr=nClust[1], ncc=nClust[2], scale="none", na.rm=F, margins=margins, main=main, xlab=NULL, ylab=NULL, zlm=limit,cexCol=2, high=colHM[1], low=colHM[2], mid=colHM[3])
                dev.off()
            }
        }
    }
    if (!is.null(samCol)) {
        for (varId in 1:length(varListAll)) {
            if (varList[varId]%in%c("os") | length(grep("t_|score",varList[varId]))==1) {
                if (outFormat=="png") {
                    png(paste("heatmapSampleColorBarLegend_",varListAll[varId],".png",sep=""),width=480,height=140)
                } else {
                    pdf(paste("heatmapSampleColorBarLegend_",varListAll[varId],".pdf",sep=""))
                }
                x=round(annColAll[,varListAll[varId]])
                if (length(grep("t_|score",varList[varId]))==1) {
                    lim=limSc
                } else {
                    lim=range(x,na.rm=T)
                    #lim=quantile(x,probs=c(.1,.9),na.rm=T)
                }
                grpUniq=lim[1]:lim[2]
                samColUniq=gray(0:(length(grpUniq)-1)/length(grpUniq))
                heatmapColorBar(limit=lim,cols=c(samColUniq[c(length(samColUniq),1)]),main=varNameAll[varId])
            } else {
                if (outFormat=="png") {
                    png(paste("heatmapSampleColorBarLegend_",varListAll[varId],".png",sep=""))
                } else {
                    pdf(paste("heatmapSampleColorBarLegend_",varListAll[varId],".pdf",sep=""))
                }
                if (varList[varId]%in%c("time")) {
                    x=annColAll[,varListAll[varId]]
                } else {
                    x=as.character(annColAll[,varListAll[varId]]); x[x==""]=NA
                }
                grpUniq=table(x)
                #grpUniq=paste(names(grpUniq)," (",grpUniq,")",sep="")
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
    if (!is.null(cloneCol)) {
        for (varId in 1:length(varFListAll)) {
            if (varFList[varId]%in%c("weight","log2FC")) {
                if (outFormat=="png") {
                    png(paste("heatmapFeatureColorBarLegend_",varFListAll[varId],".png",sep=""),width=480,height=140)
                } else {
                    pdf(paste("heatmapFeatureColorBarLegend_",varFListAll[varId],".pdf",sep=""))
                }
                x=round(annRowAll[,varFListAll[varId]])
                lim=max(abs(x),na.rm=T); lim=c(-lim,lim)
                if (varFList[varId]%in%c("log2FC")) {
                    lim=limFCmmu
                }
                grpUniq=lim[1]:lim[2]
                cloneColUniq=gray(0:(length(grpUniq)-1)/length(grpUniq))
                heatmapColorBar(limit=lim,cols=c(cloneColUniq[c(length(cloneColUniq),1)]))
            } else {
                if (outFormat=="png") {
                    png(paste("heatmapFeatureColorBarLegend_",varFListAll[varId],".png",sep=""))
                } else {
                    pdf(paste("heatmapFeatureColorBarLegend_",varFListAll[varId],".pdf",sep=""))
                }
                if (varFList[varId]%in%c("time")) {
                    x=annRowAll[,varFListAll[varId]]
                } else {
                    x=as.character(annRowAll[,varFListAll[varId]]); x[x==""]=NA
                }
                grpUniq=table(x)
                grpUniq=names(grpUniq)
                k=1:length(grpUniq)
                if (length(grpUniq)<=length(colList2)) {
                    sampleColorLegend(tls=grpUniq[k],col=colList2,legendTitle=varFNameAll[varId])
                } else if (length(grpUniq)<=length(colList)) {
                    sampleColorLegend(tls=grpUniq[k],col=colList,legendTitle=varFNameAll[varId])
                } else {
                    sampleColorLegend(tls=grpUniq[k],col=rainbow(length(grpUniq)),legendTitle=varFNameAll[varId])
                }
            }
            dev.off()
        }
    }
    
    if (is.na(nClust[1])) {
        tbl=cbind(annRow,order=1:nrow(annRow))
        write.table(tbl, paste("clusterInfoFeature",fNameOut,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)
    } else {
        if (F) {
            #png(paste("clusterVariables",fNameOut,".png",sep=""))
            pdf(paste("clusterVariables",fNameOut,".pdf",sep=""))
            plot(clustR,main=paste("Variable clusters with ",nClust[1]," main clusters marked in red",sep=""),xlab="",sub="",ylab=NULL,axes=F, cex=.2); rect.hclust(clustR,k=nClust[1])
            dev.off()
        }
        
        clustId=cutree(clustR,k=nClust[1])[clustR$order]
        k1=which(!duplicated(clustId))
        for (k in 1:length(k1)) {
            clustId[which(clustId==clustId[k1[k]])]=paste("cluster",k,sep="")
        }
        
        #tbl=as.data.frame(as.matrix(arrayData[clustR$order,]),stringsAsFactors=F)
        #tbl=data.frame(variable=clustR$labels[clustR$order],clustId,order=1:nrow(arrayData),stringsAsFactors=F)
        tbl=cbind(annRow[clustR$order,],clustId,order=1:nrow(annRow))
        write.table(tbl, paste("clusterInfoFeature",fNameOut,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)
    }
    
    if (is.na(nClust[2])) {
        tbl=cbind(annCol,order=1:nrow(annCol))
        write.table(tbl, paste("clusterInfoSample",fNameOut,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)
    } else {
        if (F) {
            #png(paste("clusterSamples",fNameOut,".png",sep=""))
            pdf(paste("clusterSamples",fNameOut,".pdf",sep=""))
            plot(clustC,main=paste("Sample clusters with ",nClust[2]," main clusters marked in red",sep=""),xlab="",sub="",ylab=NULL,axes=F, cex=.2); rect.hclust(clustC,k=nClust[2])
            dev.off()
        }
        
        clustId=cutree(clustC,k=nClust[2])[clustC$order]
        k1=which(!duplicated(clustId))
        for (k in 1:length(k1)) {
            clustId[which(clustId==clustId[k1[k]])]=paste("cluster",k,sep="")
        }
        
        tbl=cbind(annCol[clustC$order,which(!names(annCol)%in%c("order"))],clustId,order=1:nrow(annCol))
        write.table(tbl, paste("clusterInfoSample",fNameOut,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)
    }
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
if (F) {
    png("heatmaplog2FoldChangeColorBarRange.png",width=480,height=140)
    grpUniq=limFCmmu[1]:limFCmmu[2]
    colColUniq=gray(0:(length(grpUniq)-1)/length(grpUniq))
    heatmapColorBar(limit=limFCmmu,cols=c(colColUniq[c(length(colColUniq),1)],median(colColUniq)),main="log2FC")
    dev.off()
}

###########################################################
###########################################################
