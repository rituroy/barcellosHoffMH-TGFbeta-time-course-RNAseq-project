###########################################################
###########################################################

library(marray)
#source(paste(dirSrc,"functions/heatmap.5.2.R",sep=""))
#source(paste(dirSrc,"functions/heatmapAcgh.7.1.R",sep=""))
source(paste(dirSrc,"functions/heatmap.5.6.R",sep=""))
source(paste(dirSrc,"functions/heatmapAcgh.7.2.R",sep=""))

datadirG=""

outFormat="pdf"
outFormat="png"

sampleBar=""
sampleBar="cluster"

geneBar="clusterPr"
geneBar=""

centrFlag="_noCentering"
centrFlag=""

subsetFlag=""

numPr=500
pThres=10^-8
pThres=10^-6
pThres=0.05

limSc=c(-100,100)

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

colGeneId="geneId"; colIdPV="FDR"; colNamePV="QV"

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
            limFCmmu=c(-6,6)
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
                if (cohortFlag%in%c("_GSE78220")) {
                    j=which(!duplicated(phenV$patientId))
                    if (subsetFlag=="") {
                        subsetFlag=subsetName=""
                    } else {
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
                    header=paste(sub("_","",cohortFlag),subsetFlag,sep="")
                    fNameOut=paste(compFlag,cohortFlag,subsetFlag,sep="")
                    x=datV; x[datV==0]=NA
                    x=min(c(x),na.rm=T)/10
                    tbl=read.table(paste(datadirG,"geneList",cohortFlag,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
                    tbl=tbl[which(tbl$geneList==geneFlag),]
                    tbl=tbl[order(tbl$weight),]
                    iC=match(tbl$geneSym,annV$geneSym)
                    arrayData=log2(datV[iC,j]+x)
                    annRow=cbind(annV[iC,],weight=tbl$weight)
                    annCol=phenV[j,]
                    annColAll=phenV
                    varList=c("antiPd1Resp")
                    varList=c("antiPd1Resp","studySite","gender","diseaseStatus","vitalStatus","previousMapki","mutation","biopsyTime")
                    varName=paste(varList," ",sep="")
                    tbl=read.table(paste(datadirG,"scoreMat",cohortFlag,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
                    #tbl=tbl[match(phenV$id,tbl$id),which(substr(names(tbl),1,nchar("t_"))=="t_")]
                    nm1=paste("t_",geneFlag,sep="")
                    tbl=tbl[match(phenV$id,tbl$id),which(names(tbl)==nm1)]
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
                } else if (cohortFlag%in%c("_GSE33331")) {
                    header=paste(sub("_","",cohortFlag),subsetFlag,sep="")
                    fNameOut=paste(compFlag,cohortFlag,subsetFlag,sep="")
                    tbl=read.table(paste(datadirG,"geneList",cohortFlag,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
                    tbl=tbl[which(tbl$geneList==geneFlag),]
                    tbl=tbl[order(tbl$weight),]
                    iC=match(tbl$geneSym,annV$geneSym)
                    arrayData=datV[iC,]
                    annRow=cbind(annV[iC,],weight=tbl$weight)
                    annCol=phenV[,]
                    annColAll=phenV
                    varList=c("os")
                    varName=paste(varList," ",sep="")
                    if (F) {
                        tbl=read.table(paste(datadirG,"scoreMat",cohortFlag,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
                        #tbl=tbl[match(phenV$id,tbl$id),which(substr(names(tbl),1,nchar("t_"))=="t_")]
                        nm1=paste("t_",geneFlag,sep="")
                        tbl=tbl[match(phenV$id,tbl$id),which(names(tbl)==nm1)]
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
                
                k=which(varList%in%names(annCol))
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
                
                if (cohortFlag%in%c("_GSE78220","_GSE33331")) {
                    distMethod="euclidean"
                    if (cohortFlag%in%c("_GSE78220")) {
                        distMethod="pearson"
                        distMethod="cosine"
                    }
                    cloneName=rep("",nrow(annRow))
                    cloneName=annRow$geneSym
                    if (geneBar=="clusterPr") {
                        cloneCol=NULL
                    } else {
                        cloneCol=matrix(rep("white",nrow(arrayData)),nrow=1)
                        k1=1; kk=which(names(annRow)=="weight")
                        x=round(annRow[,kk]); x=x-min(x,na.rm=T)+1
                        grpUniq=sort(unique(x[!is.na(x)]))
                        limFCmmu=range(grpUniq)
                        #x=round(annRow[,kk]); x[x<limFCmmu[1]]=limFCmmu[1]; x[x>limFCmmu[2]]=limFCmmu[2]; x=x+limFCmmu[2]+1
                        grpUniq=limFCmmu[1]:limFCmmu[2]
                        cloneColUniq=gray(0:(length(grpUniq)-1)/length(grpUniq))
                        cloneCol[k1,]=cloneColUniq[x]
                        rownames(cloneCol)="weight "
                    }
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
                    samName=annCol$id2
                }
                samName=annCol$id2
                samCol=NULL
                samCol=matrix(nrow=length(varList),ncol=nrow(annCol))
                for (varId in 1:length(varList)) {
                    if (varList[varId]%in%c("os") | length(grep("t_",varList[varId]))==1) {
                        j=match(annCol$id,annColAll$id)
                        x=round(annColAll[,varList[varId]])
                        if (length(grep("t_",varList[varId]))==1) {
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
                    ncc=NA
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
                    ncr=NA
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
            if (varList[varId]%in%c("os") | length(grep("t_",varList[varId]))==1) {
                x=round(annColAll[,varListAll[varId]])
                if (length(grep("t_",varList[varId]))==1) {
                    lim=limSc
                } else {
                    lim=range(x,na.rm=T)
                    #lim=quantile(x,probs=c(.1,.9),na.rm=T)
                }
                grpUniq=lim[1]:lim[2]
                samColUniq=gray(0:(length(grpUniq)-1)/length(grpUniq))
                heatmapColorBar(limit=lim,cols=c(samColUniq[c(length(samColUniq),1)],median(samColUniq)))
            } else {
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
    if (F) {
        for (varId in 1:length(varListAll)) {
            if (outFormat=="png") {
                png(paste("heatmapSampleColorBarLegend_",varListAll[varId],".png",sep=""))
            } else {
                pdf(paste("heatmapSampleColorBarLegend_",varListAll[varId],".pdf",sep=""))
            }
            if (varList[varId]%in%c("os") | length(grep("t_",varList[varId]))==1) {
                x=round(annColAll[,varListAll[varId]])
                if (length(grep("t_",varList[varId]))==1) {
                    lim=limSc
                } else {
                    lim=range(x,na.rm=T)
                    #lim=quantile(x,probs=c(.1,.9),na.rm=T)
                }
                grpUniq=lim[1]:lim[2]
                samColUniq=gray(0:(length(grpUniq)-1)/length(grpUniq))
                heatmapColorBar(limit=lim,cols=c(samColUniq[c(length(samColUniq),1)],median(samColUniq)))
            } else {
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