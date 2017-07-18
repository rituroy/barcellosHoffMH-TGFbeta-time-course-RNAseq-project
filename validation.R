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
names(ann)[match(c("hgnc_symbol"),names(ann))]=c("geneSym")
ann$chr=as.integer(ann$chromosome_name)
ann$chr[which(ann$chromosome_name=="X")]=23
ann$chr[which(ann$chromosome_name=="Y")]=24
i=which(!is.na(ann$chr))
x=paste("chr",ann$chromosome_name,":",ann$start_position,"-",ann$end_position,sep="")[i]
write.table(x, file="ann_forUCSCliftover_hgGRCh38.txt",col.names=F,row.names=F, sep="\t",quote=F)

datadir=""
annC=read.table(paste(datadir,"ann_fromUCSCliftover_hgGRCh38toHgGRCh37.txt",sep=""),sep="\t",h=F,quote="",comment.char="",as.is=T,fill=T)
annU=read.table(paste(datadir,"ann_fromUCSCliftover_unconverted_hgGRCh38toHgGRCh37.txt",sep=""),sep="\t",h=F,quote="",comment.char="",as.is=T,fill=T)
i2=which(substr(annU[,1],1,nchar("chr"))=="chr")
x2=annU[i2,1]
x1=rep(NA,length(x))
x1[!x%in%x2]=annC[,1]
out=t(sapply(x1,function(x) {
    y=strsplit(x,":")[[1]]
    chr=y[1]
    if (!is.na(chr)) {
        if (chr=="chrX") {chr="chr23"} else if (chr=="chrY") {chr="chr24"}
    }
    chr=as.integer(sub("chr","",chr))
    pos=strsplit(y[2],"-")[[1]]
    c(chr,pos[1],pos[2])
},USE.NAMES=F))
tmp=rep("",nrow(ann))
ann2=matrix("",nrow=nrow(ann),ncol=3)
colnames(ann2)=c("chr","start","end")
ann2[i,]=out
ann2=as.data.frame(ann2,stringsAsFactors=F)
ann2$chr=sub("chr","",ann2$chr)
for (k in which(names(ann2)%in%c("chr","start","end"))) {
    ann2[,k]=as.integer(ann2[,k])
}
ann2=cbind(ann[,c("geneId","geneSym")],ann2)

annA=ann
ann2A=ann2

## -------------------
datadir="results/comparison/countNorm/donorFixedEffect/"
datadir="results/final/tables/"

stat1_4_10=read.table(paste(datadir,"stat_TGFbetaVuntreated_4hrs.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
stat1_8_10=read.table(paste(datadir,"stat_TGFbetaVuntreated_8hrs.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
stat1_12_10=read.table(paste(datadir,"stat_TGFbetaVuntreated_12hrs.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
stat1_36_10=read.table(paste(datadir,"stat_TGFbetaVuntreated_36hrs.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)

datadir="/Users/royr/Downloads/tmp/"

stat2_4_10=read.table(paste(datadir,"stat_TGFbetaVuntreated_4hrs_voom_minCnt10.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
stat2_8_10=read.table(paste(datadir,"stat_TGFbetaVuntreated_8hrs_voom_minCnt10.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
stat2_12_10=read.table(paste(datadir,"stat_TGFbetaVuntreated_12hrs_voom_minCnt10.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
stat2_36_10=read.table(paste(datadir,"stat_TGFbetaVuntreated_36hrs_voom_minCnt10.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)


stat1=stat1_12_10; stat2=stat2_12_10
stat1=stat1_8_10; stat2=stat2_8_10
stat1=stat1_4_10; stat2=stat2_4_10
stat1=stat1_36_10; stat2=stat2_36_10
k=match(names(stat1),names(stat2)); k1=which(!is.na(k)); k2=k[k1]
for (k in 1:length(k1)) {
    if (any(is.na(stat1[,k1[k]])!=is.na(stat2[,k2[k]])) | any(stat1[,k1[k]]!=stat2[,k2[k]],na.rm=T)) {
        cat(k,names(stat1)[k1[k]],"\n")
        print(max(abs(stat1$logCPM-stat2$logCPM),na.rm=T)/min(abs(stat1$logCPM),na.rm=T))
    }
}

stat1_4_10$t=stat2_4_10$t
stat1_8_10$t=stat2_8_10$t
stat1_12_10$t=stat2_12_10$t
stat1_36_10$t=stat2_36_10$t

## -------------------
## Weight

#library(pls)

ann=annA
ann2=ann2A

colId=c("geneId","logFC","logCPM","PValue","Qvalue")
colId=c("geneId","logFC","logCPM","t","PValue","Qvalue")
pThres=0.05
geneList=c("TGFbetaVuntreated_8hrs_qv0.05","TGFbetaVuntreated_8hrs8fold_qv0.05","TGFbetaVuntreated_8hrs16fold_qv0.05","TGFbetaVuntreated_8hrs_12hrs_qv0.05_sameDir","TGFbetaVuntreated_8hrs2fold_12hrs_qv0.05_sameDir","TGFbetaVuntreated_8hrs2fold_12hrs2fold_qv0.05_sameDir","TGFbetaVuntreated_8hrs8fold_12hrs8fold_qv0.05_sameDir","TGFbetaVuntreated_12hrs_qv0.05","TGFbetaVuntreated_12hrs8fold_qv0.05","TGFbetaVuntreated_36hrs_qv0.05","TGFbetaVuntreated_36hrs8fold_qv0.05")
for (geneFlag in geneList) {
    cat("\n\n============ ",geneFlag," ============\n",sep="")
    switch(geneFlag,
    "TGFbetaVuntreated_8hrs_qv0.05"={
        lfcThres=0
        stat2=stat1_8_10[,colId]
        i1=which(stat2$Qvalue<pThres)
    },
    "TGFbetaVuntreated_8hrs8fold_qv0.05"={
        lfcThres=3
        stat2=stat1_8_10[,colId]
        i1=which(stat2$Qvalue<pThres & abs(stat2$logFC)>=lfcThres)
    },
    "TGFbetaVuntreated_8hrs16fold_qv0.05"={
        lfcThres=4
        stat2=stat1_8_10[,colId]
        i1=which(stat2$Qvalue<pThres & abs(stat2$logFC)>=lfcThres)
    },
    "TGFbetaVuntreated_8hrs_12hrs_qv0.05_sameDir"={
        lfcThres=0
        i=match(stat1_8_10$geneId,stat1_12_10$geneId); i1=which(!is.na(i)); i2=i[i1]
        stat2=stat1_8_10[i1,colId]
        stat2=stat2[which(stat1_12_10$Qvalue[i2]<pThres & sign(stat2$logFC)==sign(stat1_12_10$logFC[i2])),]
        i1=which(stat2$Qvalue<pThres)
    },
    "TGFbetaVuntreated_8hrs2fold_12hrs_qv0.05_sameDir"={
        lfcThres=1
        i=match(stat1_8_10$geneId,stat1_12_10$geneId); i1=which(!is.na(i)); i2=i[i1]
        stat2=stat1_8_10[i1,colId]
        stat2=stat2[which(stat1_12_10$Qvalue[i2]<pThres & sign(stat2$logFC)==sign(stat1_12_10$logFC[i2])),]
        i1=which(stat2$Qvalue<pThres)
    },
    "TGFbetaVuntreated_8hrs2fold_12hrs2fold_qv0.05_sameDir"={
        lfcThres=1
        i=match(stat1_8_10$geneId,stat1_12_10$geneId); i1=which(!is.na(i)); i2=i[i1]
        stat2=stat1_8_10[i1,colId]
        stat2=stat2[which(stat1_12_10$Qvalue[i2]<pThres & sign(stat2$logFC)==sign(stat1_12_10$logFC[i2]) & abs(stat2$logFC)>=lfcThres & abs(stat1_12_10$logFC[i2])>=lfcThres),]
        i1=which(stat2$Qvalue<pThres)
    },
    "TGFbetaVuntreated_8hrs8fold_12hrs8fold_qv0.05_sameDir"={
        lfcThres=3
        i=match(stat1_8_10$geneId,stat1_12_10$geneId); i1=which(!is.na(i)); i2=i[i1]
        stat2=stat1_8_10[i1,colId]
        stat2=stat2[which(stat1_12_10$Qvalue[i2]<pThres & sign(stat2$logFC)==sign(stat1_12_10$logFC[i2]) & abs(stat2$logFC)>=lfcThres & abs(stat1_12_10$logFC[i2])>=lfcThres),]
        i1=which(stat2$Qvalue<pThres)
    },
    "TGFbetaVuntreated_12hrs_qv0.05"={
        lfcThres=0
        stat2=stat1_12_10[,colId]
        i1=which(stat2$Qvalue<pThres)
    },
    "TGFbetaVuntreated_12hrs8fold_qv0.05"={
        lfcThres=3
        stat2=stat1_12_10[,colId]
        i1=which(stat2$Qvalue<pThres & abs(stat2$logFC)>=lfcThres)
    },
    "TGFbetaVuntreated_36hrs_qv0.05"={
        lfcThres=0
        stat2=stat1_36_10[,colId]
        i1=which(stat2$Qvalue<pThres)
    },
    "TGFbetaVuntreated_36hrs8fold_qv0.05"={
        lfcThres=3
        stat2=stat1_36_10[,colId]
        i1=which(stat2$Qvalue<pThres & abs(stat2$logFC)>=lfcThres)
    }
    )
    stat2=cbind(stat2,ann2[match(stat2$geneId,ann2$geneId),which(!names(ann2)%in%names(stat2))])
    cat("No. of genes: ",length(i1),"\n",sep="")
    wtMat=matrix(nrow=length(i1), ncol=3,dimnames=list(stat2$geneId[i1],c("t","pca","plsr")))
    wtMat[,"t"]=stat2$t[i1]
    if (F) {
        cat("===============",colnames(log2FoldChng)[datasetId],"\n")
        resp=rep(NA,length(grpName))
        resp[which(grpName%in%strsplit(grpName1[grpId],".",fixed=T)[[1]])]=0
        resp[which(grpName%in%strsplit(grpName2[grpId],".",fixed=T)[[1]])]=1
        resp.na=!is.na(resp)
        resp=resp[resp.na]
        print(table(resp))
        x=apply(expr[clId1,],1,function(x) {x-median(x,na.rm=T)})
        colnames(x)=paste("X",1:ncol(x),sep="")
        thisData=as.data.frame(cbind(resp,x[resp.na,]))
        
        if (min(table(resp))>1) {
            fmla=as.formula(paste(" ~ ", paste(colnames(x), collapse= "+")))
            fit=prcomp(fmla, center=F, scale=F, data = thisData)
            wtMat$pca[,datasetId]=fit$rotation[,1]
            
            fit=plsr(resp ~ ., ncomp=1, method="kernelpls", center=F, scale=F, data = thisData)
            wtMat$plsr[,datasetId]=coef(fit,comp=1)[,1,1]
            
        }
    }
    tbl=cbind(geneId=stat2$geneId[i1],as.data.frame(wtMat))
    write.table(tbl, paste("wtMat_",geneFlag,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)
}

#####################################################################################
#####################################################################################
## Sample score

getScore=function(logFlag=F,centerFlag=T,scaleFlag=F) {
    signatFlag=c("_uniqGeneSym","_scale")
    ann=annA
    
    sdVec=apply(datV,1,sd,na.rm=T)
    i2=order(sdVec,decreasing=T)
    datadir=""
    wtList=c("t","pca","plsr")
    scoreMat=matrix(nrow=ncol(datV),ncol=3*length(geneList),dimnames=list(phenV$id,paste(wtList,"_",rep(geneList,each=length(wtList)),sep="")))
    tmp=rep(NA,length(geneList))
    tmpC=geneList
    geneInfo=data.frame(geneList,train=tmp,val=tmp,stringsAsFactors=F)
    geneAnno=NULL
    for (geneFlag in geneList) {
        wtMat=read.table(paste(datadir,"wtMat_",geneFlag,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
        x=as.integer(strsplit(strsplit(geneFlag,"_")[[1]][2],"hrs")[[1]][1])
        if (x==8) {
            stat2=stat1_8_10
        } else if (x==12) {
            stat2=stat1_12_10
        } else if (x==36) {
            stat2=stat1_36_10
        } else {
            stat2=NULL
        }
        statThis=stat2[match(wtMat$geneId,stat2$geneId),]
        i1=match(wtMat$geneId,ann$geneId)
        i=match(tolower(ann$geneSym[i1]),tolower(annV$geneSym[i2])); i11=which(!is.na(i)); i21=i2[i[!is.na(i)]]
        length(i21)
        
        k=which(geneInfo$geneList==geneFlag)
        geneInfo$train[k]=length(i1)
        geneInfo$val[k]=length(i21)
        geneAnno=rbind(geneAnno,matrix(c(rep(geneFlag,length(i21)),annV$geneSym[i21],wtMat$t[i11],statThis$logFC[i11]),ncol=4,byrow=F))
        
        x=datV[i21,]
        if (logFlag) {
            x0=x
            x[x0==0]=NA
            x=min(c(x),na.rm=T)/10
            x=log2(x0+x)
        }
        if (centerFlag) {
                for (i in 1:nrow(x)) {
                x[i,]=x[i,]-median(x[i,],na.rm=T)
            }
        }
        if (scaleFlag) {
            #centr=apply(x,1,mad,na.rm=T)
            for (i in 1:nrow(x)) {
                x[i,]=x[i,]/mad(x[i,],na.rm=T)
            }
            #cat("\n\n=============== ",geneFlag," ==================\n")
            #print(table(centr==0))
        }
        score=rep(0,ncol(x))
        for (j in 1:ncol(x)) {
            y=wtMat$t[i11]*x[,j]
            ii=which(!is.infinite(y))
            score[j]=sum(y[ii],na.rm=T)
        }
        if (signatFlag[2]=="_scale") score=100*score/max(abs(score))
        scoreMat[,which(colnames(scoreMat)==paste("t_",geneFlag,sep=""))]=score
    }
    tbl=as.data.frame(scoreMat)
    tbl=cbind(id=phenV$id,tbl)
    write.table(tbl, paste("scoreMat",cohortFlag,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)

    colnames(geneAnno)=c("geneList","geneSym","weight","logFC")
    write.table(geneAnno, paste("geneList",cohortFlag,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)

    write.table(geneInfo, paste("geneInfo",cohortFlag,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)
    
    invisible(scoreMat)
}

#####################################################################################
#####################################################################################
cohortFlag="_GSE33331"

## -------------------
datadir="docs/GSE33331/"
fName="GSE33331_series_matrix.txt"
tmp=read.table(paste(datadir,fName,sep=""),sep="\n",h=F,quote="",comment.char="",as.is=T,fill=T,nrow=200)
k=which(substr(tmp[,1],1,nchar("!Sample_title"))=="!Sample_title")
id=gsub("\"","",strsplit(unlist(tmp[k,]),"\t")[[1]][-1])
k=which(substr(tmp[,1],1,nchar("!Sample_characteristics_ch1"))=="!Sample_characteristics_ch1")
phenV=read.table(paste(datadir,fName,sep=""),sep="\t",h=F,quote="",comment.char="",as.is=T,fill=T,skip=k[1]-1,nrow=2*length(k))
phenV=apply(t(as.matrix(phenV[which(phenV[,1]=="!Sample_characteristics_ch1"),-1])),2,function(x) {
    y=gsub("\"","",x)
    nm=y[which(y!="")]
    nm=strsplit(nm[1],": ")[[1]][1]
    z=sapply(y,function(y) strsplit(y,": ")[[1]][2],USE.NAMES=F)
    c(nm,z)
})
nm=phenV[1,]
colnames(phenV)[match(c("tissue","genotype","survival (months)"),nm)]=c("tissue","genotype","os")
rownames(phenV)=NULL
phenV=as.data.frame(cbind(id=id,phenV[-1,]),stringsAsFactors=F)
k=which(substr(tmp[,1],1,nchar("!series_matrix_table_begin"))=="!series_matrix_table_begin")
k=68
datV=read.table(paste(datadir,fName,sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,skip=k,nrow=-1)
datadir="/Users/royr/Downloads/HG-U133_Plus_2-na36-annot-csv/"
fName="HG-U133_Plus_2.na36.annot.csv"
if (F) {
    tmp=read.table(paste(datadir,fName,sep=""),sep="\n",h=F,quote="",comment.char="",as.is=T,fill=T,nrow=100)
    k=which(substr(tmp[,1],1,nchar("\"Probe Set ID\""))=="\"Probe Set ID\"")
    id=gsub("\"","",strsplit(unlist(tmp[k,]),",")[[1]])
    if (T) {
        annV=read.table(paste(datadir,fName,sep=""),sep=",",h=F,quote="",comment.char="",as.is=T,fill=T,skip=26,nrow=-1)
        names(annV)[match(c("Probe Set ID","GeneChip Array","Species Scientific Name","Annotation Date","Sequence Type","Sequence Source","Transcript ID(Array Design)","Target Description","Representative Public ID","Archival UniGene Cluster","UniGene ID","Genome Version","Alignments","Gene Title","Gene Symbol","Chromosomal Location","Unigene Cluster Type","Ensembl","Entrez Gene","SwissProt","EC","OMIM","RefSeq Protein ID","RefSeq Transcript ID","FlyBase","AGI","WormBase","MGI Name","RGD Name","SGD accession number","Gene Ontology Biological Process","Gene Ontology Cellular Component","Gene Ontology Molecular Function","Pathway","InterPro","Trans Membrane","QTL","Annotation Description","Annotation Transcript Cluster","Transcript Assignments","Annotation Notes"),id)]=
            c("id","GeneChip Array","Species Scientific Name","Annotation Date","Sequence Type","Sequence Source","Transcript ID(Array Design)","Target Description","Representative Public ID","Archival UniGene Cluster","UniGene ID","Genome Version","Alignments","Gene Title","geneSym","Chromosomal Location","Unigene Cluster Type","Ensembl","Entrez Gene","SwissProt","EC","OMIM","RefSeq Protein ID","RefSeq Transcript ID","FlyBase","AGI","WormBase","MGI Name","RGD Name","SGD accession number","Gene Ontology Biological Process","Gene Ontology Cellular Component","Gene Ontology Molecular Function","Pathway","InterPro","Trans Membrane","QTL","Annotation Description","Annotation Transcript Cluster","Transcript Assignments","Annotation Notes")
        for (k in 1:ncol(annV)) {
            annV[,k]=gsub("\"","",annV[,k])
        }
        i=order(annV$id)
        cbind(i,annV$id[i])[1:5,]
        annV=annV[,c("id","geneSym")]
    }
    if (F) {
        annV=read.table(paste(datadir,fName,sep=""),sep="\n",h=F,quote="",comment.char="",as.is=T,fill=T,skip=26,nrow=10)
        annV=t(sapply(annV[,1],function(x) {
            y=strsplit(x,",")[[1]][c(1,15)]
            y=gsub("\"","",y)
        },USE.NAMES=F))
    }
}
fName="HG-U133_Plus_2.na36.annot_RRedit_lean.csv"
annV=read.table(paste(datadir,fName,sep=""),sep=",",h=T,quote="",comment.char="",as.is=T,fill=T)
names(annV)[match(c("Probe.Set.ID","Gene.Symbol","Chromosomal.Location","Unigene.Cluster.Type","Ensembl","Entrez.Gene"),names(annV))]=c("id","geneSym","chrLoc","unigeneClustType","ensemblGeneId","entrezGeneId")

phenV$geoAcc=sapply(names(datV)[-1],function(x) {strsplit(x,".",fixed=T)[[1]][2]},USE.NAMES=F)
for (k in which(names(phenV)%in%c("os"))) {
    phenV[,k]=as.numeric(phenV[,k])
}
id=gsub("\"","",datV[,1])
id[!id%in%annV$id]
table(tolower(annV$geneSym)%in%tolower(annA$geneSym))
i=match(id,annV$id); i1=which(!is.na(i)); i2=i[i1]
datV=as.matrix(datV[i1,-1])
colnames(datV)=phenV$id
rownames(datV)=id[i1]
annV=annV[i2,]

y=c()
for (k in 1:ncol(phenV)) {
    x=sum(!duplicated(phenV[!is.na(phenV[,k]),k]))
    if (x>1 & x<=10) {
        cat("============= ",k,": ",names(phenV)[k],"\n")
        print(table(phenV[,k]))
        y=c(y,names(phenV)[k])
    }
}
paste(y,collapse=",")

#scoreMat=getScore()
#scoreMat=getScore(logFlag=F,centerFlag=T,scaleFlag=T)
#scoreMat=getScore(logFlag=F,centerFlag=T,scaleFlag=F)
scoreMat=getScore(logFlag=F,centerFlag=F,scaleFlag=F)

datObj31=list(dat=datV,phen=phenV,ann=annV,score=scoreMat)

## -------------------

## -------------------
cohortFlag="_GSE78220"

datadir="docs/GSE78220/"

fName="GSE78220_series_matrix.txt"
fName="GSE78220_series_matrix_RRedit.txt"
tmp=read.table(paste(datadir,fName,sep=""),sep="\n",h=F,quote="",comment.char="",as.is=T,fill=T)
k=which(substr(tmp[,1],1,nchar("!Sample_title"))=="!Sample_title")
id=strsplit(unlist(tmp[k,]),"\t")[[1]][-1]
phenV=read.table(paste(datadir,fName,sep=""),sep="\t",h=F,quote="",comment.char="",as.is=T,fill=T,skip=k[1]-1,nrow=2*length(k))
k=which(substr(tmp[,1],1,nchar("!Sample_characteristics_ch1"))=="!Sample_characteristics_ch1")
phenV=read.table(paste(datadir,fName,sep=""),sep="\t",h=F,quote="",comment.char="",as.is=T,fill=T,skip=k[1]-1,nrow=2*length(k))
phenV=apply(t(as.matrix(phenV[which(phenV[,1]=="!Sample_characteristics_ch1"),-1])),2,function(x) {
    y=gsub("\"","",x)
    nm=y[which(y!="")]
    nm=strsplit(nm[1],": ")[[1]][1]
    if (nm%in%c("braf","nras","nf1")) {
        nm="mutation"
        y=paste(nm,": ",y,sep="")
    }
    z=sapply(y,function(y) strsplit(y,": ")[[1]][2],USE.NAMES=F)
    c(nm,z)
})
nm=phenV[1,]
colnames(phenV)=phenV[1,]
colnames(phenV)[match(c("patient id","anti-pd-1 response","study site","gender","age (yrs)","disease status","overall survival (days)","vital status","previous mapki","anatomical location","mutation","treatment","biopsy time","tissue"),colnames(phenV))]=
c("patientId","antiPd1Resp","studySite","gender","age","diseaseStatus","overallSurvival","vitalStatus","previousMapki","anatomicalLocation","mutation","treatment","biopsyTime","tissue")
rownames(phenV)=NULL
phenV=as.data.frame(cbind(id=id,phenV[-1,]),stringsAsFactors=F)
datV=read.table(paste(datadir,"GSE78220_PatientFPKM.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,nrow=-1)
annV=data.frame(id=sapply(datV$Gene,function(x) {gsub("-","_",x)},USE.NAMES=F),geneSym=datV$Gene,stringsAsFactors=F)
id=sapply(names(datV),function(x) {strsplit(x,".",fixed=T)[[1]][1]},USE.NAMES=F)
datV=as.matrix(datV[,match(phenV$id,id)])
colnames(datV)=phenV$id

y=c()
for (k in 1:ncol(phenV)) {
    x=sum(!duplicated(phenV[!is.na(phenV[,k]),k]))
    if (x>1 & x<=10) {
        cat("============= ",k,": ",names(phenV)[k],"\n")
        print(table(phenV[,k]))
        y=c(y,names(phenV)[k])
    }
}
paste(y,collapse=",")

#scoreMat=getScore(logFlag=T,centerFlag=T,scaleFlag=T)
#scoreMat=getScore(logFlag=T,centerFlag=T,scaleFlag=F)
scoreMat=getScore(logFlag=T,centerFlag=F,scaleFlag=F)

datObj20=list(dat=datV,phen=phenV,ann=annV,score=scoreMat)

#####################################################################################
#####################################################################################
## DE genes
## Nothing significant for _GSE78220, antiPd1RespBi

library(limma)
library(qvalue)

cohortFlag="_GSE78220"
subsetFlag="_biopsyPreTreat"

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
}
phenThis=cbind(phenV[j,],scoreMat[j,grep("t_",colnames(scoreMat))])
phenThis$antiPd1RespBi=phenThis$antiPd1Resp
phenThis$antiPd1RespBi[which(phenThis$antiPd1Resp%in%c("Complete Response","Partial Response"))]="CR/PR"
i=apply(datV,1,function(x) mean(x,na.rm=T)!=0)
i=apply(datV,1,function(x) mean(x,na.rm=T)!=0 & (sum(x==0,na.rm=T)/sum(!is.na(x)))<0.05)
dat=log2(datV[i,j]+1)
ann=annV[i,]

varList="antiPd1RespBi"
grp=phenThis[,varList]
j=which(!is.na(grp))
design=model.matrix(~grp[j])
fit=lmFit(dat[,j],design)
fit=eBayes(fit)
colId=2
stat2=ann
names(stat2)=c("geneId","geneSym")
stat2=cbind(stat2,logFC=fit$coef[,colId],t=fit$t[,colId],pv=fit$p.value[,colId])
stat2$qv=NA
i=which(!is.na(stat2$pv))
stat2$qv[i]=p.adjust(stat2$pv[i],method="holm")
stat20_1=stat2

pThres=0.05
png(paste("plots",cohortFlag,".png",sep=""))
par(mfrow=c(2,2))
hist(stat2$pv)
i=which(stat2$qv<pThres)
plot(stat2$logFC,-log10(stat2$pv),main=paste(varList,"\nNo. with qv<",pThres," = ",sum(i),sep=""),pch=19,cex=.8)
points(stat2$logFC[i],-log10(stat2$pv[i]),pch=19,cex=.8,col="red")
qqnorm(stat2$pv); qqline(stat2$pv)
dev.off()
dim(dat)

#####################################################################################
#####################################################################################
## Association of sample score with clinical variables

cohortFlag="_GSE78220"
cohortFlag="_GSE33331"

pThres=0.05

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

tbl=NULL
if (cohortFlag=="_GSE33331") {
    library(survival)

    subsetFlag=""
    phenThis=cbind(phenV,event=rep(1,nrow(phenV)),scoreMat[,grep("t_",colnames(scoreMat))])
    for (k in grep("t_",names(phenThis))) {
        #cat("\n\n=========== ",names(phenThis)[k]," ===========\n")
        x=phenThis[,k]
        res=coxph(Surv(os,event)~x,data=phenThis)
        #print(summary(res)$coef)
        tbl2=c(subsetFlag,"os",names(phenThis)[k],summary(res)$coef[1,"Pr(>|z|)"],"coxph")
        tbl=rbind(tbl,tbl2)
    }
    rownames(tbl)=NULL
    tbl=as.data.frame(tbl,stringsAsFactors=F)
    names(tbl)=c("subset","variable","scoreType","pv","testType")
    tbl$pv=as.numeric(tbl$pv)
}


## -------------------
if (cohortFlag=="_GSE78220") {
    
    print(table(phenV$antiPd1Resp,phenV$studySite))
    print(table(phenV$antiPd1Resp,phenV$biopsyTime))
    j=which(phenV$biopsyTime=="pre-treatment")
    x=table(phenV$antiPd1Resp[j],phenV$studySite[j])
    print(x)
    cat("Fisher's exact test pv ",fisher.test(x)$p.value,"\n",sep="")
    
    library(coin)

    subsetFlag="_ucla"
    subsetFlag=""
    subsetFlag="_biopsyPreTreat"
    subsetFlag="_ucla_biopsyPreTreat"
    subsetList=c("","_ucla","_biopsyPreTreat","_ucla_biopsyPreTreat")
    for (subsetFlag in subsetList) {
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
        }

        phenThis=cbind(phenV[j,],scoreMat[j,grep("t_",colnames(scoreMat))])
        phenThis$antiPd1RespBi=phenThis$antiPd1Resp
        phenThis$antiPd1RespBi[which(phenThis$antiPd1Resp%in%c("Complete Response","Partial Response"))]="CR/PR"
        varList=c("antiPd1Resp","studySite","gender","diseaseStatus","vitalStatus","previousMapki","mutation","biopsyTime")
        varList=c("antiPd1Resp","gender","diseaseStatus","vitalStatus","previousMapki","mutation")
        varList=c("antiPd1Resp","antiPd1RespBi","gender","diseaseStatus","vitalStatus","previousMapki","mutation")
        for (k1 in match(varList,names(phenThis))) {
            #cat("\n\n==== ",names(phenThis)[k1]," ===========\n")
            for (k in grep("t_",names(phenThis))) {
                #cat("\n\n=========== ",names(phenThis)[k]," ===========\n")
                x=phenThis[,k]
                y=table(phenThis[,k1])
                if (length(y)==2) {
                    testType="wilcox"
                    res=wilcox_test(x~as.factor(phenThis[,k1]))
                } else {
                    testType="kruskal"
                    res=kruskal_test(x~as.factor(phenThis[,k1]))
                }
                #print(pvalue(res))
                pv=pvalue(res)
                tbl2=c(subsetFlag,names(phenThis)[k1],names(phenThis)[k],pv,testType)
                tbl=rbind(tbl,tbl2)
                if (pv<pThres) {
                    png(paste("boxplot_sampleScore_",names(phenThis)[k1],"_",names(phenThis)[k],cohortFlag,subsetFlag,".png",sep=""))
                    boxplot(x~as.factor(phenThis[,k1]),main=paste(cohortFlag,ifelse(subsetFlag=="","",", "),subsetFlag,"\n",sub("t_","",names(phenThis)[k]),": pv ",signif(pv,2),sep=""),ylab="Sample score")
                    dev.off()
                    if (names(phenThis)[k1]=="antiPd1RespBi") {
                        k2=which(names(phenThis)=="antiPd1Resp")
                        y=table(phenThis[,k2])
                        if (length(y)==2) {
                            res=wilcox_test(x~as.factor(phenThis[,k2]),exact=T)
                        } else {
                            res=kruskal_test(x~as.factor(phenThis[,k2]))
                        }
                        pv=pvalue(res)
                        png(paste("boxplot_sampleScore_",names(phenThis)[k2],"_",names(phenThis)[k],cohortFlag,subsetFlag,".png",sep=""))
                        boxplot(x~as.factor(phenThis[,k2]),main=paste(cohortFlag,ifelse(subsetFlag=="","",", "),subsetFlag,"\n",sub("t_","",names(phenThis)[k]),": pv ",signif(pv,2),sep=""),ylab="Sample score")
                        dev.off()
                    }
                }
            }
        }
    }
    rownames(tbl)=NULL
    tbl=as.data.frame(tbl,stringsAsFactors=F)
    names(tbl)=c("subset","variable","scoreType","pv","testType")
    tbl$pv=as.numeric(tbl$pv)
}

## -------------------
cat(cohortFlag,"\n",sep="")
tbl2=tbl
tbl2$pv=signif(tbl2$pv,2)
tbl2[tbl$pv<pThres,]
"
GSE78220

>     print(table(phenV$antiPd1Resp,phenV$studySite))

                     UCLA VIC
Complete Response      5   0
Partial Response       7   3
Progressive Disease   11   2

UCLA VIC
Complete Response      5   0
Partial Response       7   3
Progressive Disease   10   2
Fisher's exact test pv 0.5571659

----------

logFlag=T, centerFlag=F, scaleFlag=F
> tbl[tbl$pv<.05,]
subset      variable                                               scoreType     pv testType
56                  previousMapki                         t_TGFbetaVuntreated_8hrs_qv0.05 0.0034   wilcox
59                  previousMapki           t_TGFbetaVuntreated_8hrs_12hrs_qv0.05_sameDir 0.0084   wilcox
60                  previousMapki      t_TGFbetaVuntreated_8hrs2fold_12hrs_qv0.05_sameDir 0.0084   wilcox
61                  previousMapki t_TGFbetaVuntreated_8hrs2fold_12hrs2fold_qv0.05_sameDir 0.0130   wilcox
63                  previousMapki                        t_TGFbetaVuntreated_12hrs_qv0.05 0.0450   wilcox
133           _ucla previousMapki                         t_TGFbetaVuntreated_8hrs_qv0.05 0.0480   wilcox
168 _biopsyPreTreat antiPd1RespBi                   t_TGFbetaVuntreated_8hrs16fold_qv0.05 0.0400   wilcox
174 _biopsyPreTreat antiPd1RespBi                   t_TGFbetaVuntreated_12hrs8fold_qv0.05 0.0450   wilcox
210 _biopsyPreTreat previousMapki                         t_TGFbetaVuntreated_8hrs_qv0.05 0.0064   wilcox
213 _biopsyPreTreat previousMapki           t_TGFbetaVuntreated_8hrs_12hrs_qv0.05_sameDir 0.0160   wilcox
214 _biopsyPreTreat previousMapki      t_TGFbetaVuntreated_8hrs2fold_12hrs_qv0.05_sameDir 0.0160   wilcox
215 _biopsyPreTreat previousMapki t_TGFbetaVuntreated_8hrs2fold_12hrs2fold_qv0.05_sameDir 0.0240   wilcox
220 _biopsyPreTreat previousMapki                   t_TGFbetaVuntreated_36hrs8fold_qv0.05 0.0310   wilcox
"

#####################################################################################
#####################################################################################
