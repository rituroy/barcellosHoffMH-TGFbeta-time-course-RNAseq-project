R.utils::use("aroma.seq, edgeR, aroma.light, matrixStats")

#fName="/compbio/data/annotationData/organisms/Homo_sapiens/GRCh38,hg38/Ensembl/79/Homo_sapiens.GRCh38.79.gtf"

dataset <- "barcellosHoffMHTgfbTC"
dataset <- "tmp"

cohort="hg"
organism <- "Homo_sapiens"

#counts <- HTSeqCountDataSet$byName(dataset, tags = "tophat2,pe,gtf", organism = organism)
counts <- HTSeqCountDataSet$byName(dataset, tags = "tophat2,gtf", organism = organism)
counts <- setFullNamesTranslator(counts, function(name, ...) chartr(";", "/", name))
counts

data <- extractMatrix(counts, column = 2L, colClass = "integer")
dge <- extractDGEList(counts)

db <- TabularTextFile("data/sampleInfo.txt")
samples <- readDataFrame(db)
samples=samples[match(colnames(data),samples$id),]
samples1=samples

data2=data
rho <- cor(data2, method = "spearman")
round(summary(as.vector(rho)),2)
"
Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
0.2550  0.8979  0.9092  0.8159  0.9171  1.0000
"
j=which(samples$time==4)
round(rho[j,j],2)

rho <- cor(data, method = "spearman")
round(summary(as.vector(rho)),2)

dge$samples[, -1]
samples=cbind(samples,dge$samples[,c("lib.size","norm.factors")])
write.table(cbind(geneId=rownames(dge$counts),dge$counts), file=paste("count_raw_",organism,".txt",sep=""),col.names=T,row.names=F, sep="\t",quote=F)
write.table(samples, file=paste("sample_",organism,".txt",sep=""),col.names=T,row.names=F, sep="\t",quote=F)


dge$samples$group <- samples$treat

dge1=dge

#--------------------------
dge=dge1
dge <- calcNormFactors(dge, method = "TMM")
j=j[which(!rownames(dge$samples)[j]%in%c("Donor_1_TGFbeta_4hrs","Donor_3_4_hrs"))]
lcpmT=cpm(dge,log=TRUE)
tbl=as.data.frame(lcpmT)
tbl=cbind(ensembl_gene_id=rownames(tbl),tbl)
write.table(tbl, file=paste("logCPM.txt",sep=""),col.names=T,row.names=F, sep="\t",quote=F)
gzip("logCPM_allSamples.txt")

#--------------------------

dge=dge1

#dge=dge[,which(rownames(dge$samples)!="sample6")]

if (F) {
    for (samId in rownames(dge1$samples)) {
        cat("\n\n===========",samId,"=============\n\n")
        dge=dge1[,which(rownames(dge1$samples)!=samId)]

        samples=samples1[match(rownames(dge$samples),samples1$id),]

        colSam=rainbow(ncol(dge$counts))

        dge <- calcNormFactors(dge, method = "TMM")

        png("densityPlot_normCount.png")
        par(mfcol=c(3,3))
        xlim=range(c(dge$counts),na.rm=T)
        xlim=c(0,50)
        j=1
        plot(density(dge$counts[,j],na.rm=T),xlim=xlim,main=paste("Densities of scale-normalized log2 gene counts across the ",ncol(dge$counts)," samples",sep=""),xlab="log(count)",col=colSam[j])
        for (j in 2:ncol(dge$counts)) {
        #	lines(density(dge$counts[,j],na.rm=T),col=colSam[j])
            plot(density(dge$counts[,j],na.rm=T),col=colSam[j],xlim=xlim)
        }
        dev.off()

        dge$samples$group <- samples$treat

        dgeT <- dge[, !is.na(dge$samples$group)]

        dgeT <- estimateDisp(dgeT, trend = "none", robust = TRUE)
        print("table(is.infinite(dgeT$prior.df))")
        print(table(is.infinite(dgeT$prior.df)))
        "
        FALSE  TRUE
        30728 34488
        "
        
        save(dgeT,file="dgeT_tmp.RData")
    }
}

if (F) {
    ## NOT USED
    ## Use the annotation from analysis.R
    #source("biomartify.1.2.R")
    source("/home/royr/shared/biomartify.R")
    y=apply(dgeT$counts,1,function(x) {
        mean(log2(x),na.rm=T)
    })
    ann <- biomartify(data.frame(geneId=rownames(dgeT$counts),logFC=y,logCPM=y,PValue=y,FDR=y,stringsAsFactors=F), organism = ifelse(organism=="Mus_musculus","MusMusculus","HomoSapiens"))
    ann <- ann[,which(!names(ann)%in%c("logFC","logCPM","PValue","FDR"))]
    write.table(ann,paste("ann_",organism,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)
}

## -----------------------------------------
subsetFlag=""

subsetList=paste("_",sort(unique(samples1$time)),"hrs",sep="")

for (subsetFlag in subsetList) {
    fName1=paste("_",organism,subsetFlag,sep="")
    dge=dge1
    samples=samples1[match(rownames(dge$samples),samples1$id),]

    j=1:nrow(dge$samples)
    if (subsetFlag!="") j=which(paste("_",samples$time,"hrs",sep="")==subsetFlag)
    subsetName=subsetFlag
    dge=dge[,j]

    colSam=rainbow(ncol(dge$counts))

    dge <- calcNormFactors(dge, method = "TMM")
    samples=samples1[match(rownames(dge$samples),samples1$id),]
    samples=cbind(samples,dge$samples[,c("lib.size","norm.factors")])

    datadir="results/comparison/"
    datadir=""
    ann <- read.table(paste(datadir,"ann_",organism,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
    ann$length=sapply(ann$chromosome,function(x) {
                      y=strsplit(x,":")[[1]][2]
                      y=strsplit(y," ")[[1]][1]
                      y=diff(as.integer(strsplit(y,"-")[[1]]))
                      y
                      },USE.NAMES=F)
    i=match(rownames(dge$counts),sub("*","",ann$gene_id,fixed=T))
    x=rpkm(dge, gene.length=ann$length[i], normalized.lib.sizes=TRUE, log=FALSE, prior.count=0.25)
    write.table(cbind(geneId=rownames(dge$counts),x), file=paste("rpkm",fName1,".txt",sep=""),col.names=T,row.names=F, sep="\t",quote=F)
    
#	save(dge,file=paste("dge",fName1,".RData",sep=""))
    write.table(cbind(geneId=rownames(dge$counts),dge$counts), file=paste("count",fName1,".txt",sep=""),col.names=T,row.names=F, sep="\t",quote=F)
    write.table(samples, file=paste("sample",fName1,".txt",sep=""),col.names=T,row.names=F, sep="\t",quote=F)

    if (F) {
        png("densityPlot_normCount.png",width=3*240, height=3*240)
        par(mfcol=c(3,3))
        xlim=range(c(dge$counts),na.rm=T)
        xlim=c(0,50)
        xlim=NULL
        j=1
        plot(density(dge$counts[,j],na.rm=T),xlim=xlim,main=paste("Densities of scale-normalized log2 gene counts across the ",ncol(dge$counts)," samples",sep=""),xlab="log(count)",col=colSam[j])
        for (j in 2:ncol(dge$counts)) {
            lines(density(dge$counts[,j],na.rm=T),col=colSam[j])
        #	plot(density(dge$counts[,j],na.rm=T),col=colSam[j],xlim=xlim)
        }
        dev.off()
    }

    if (F) {
        png(paste("densityPlot_normCount",subsetName,".png",sep=""),width=3*240, height=3*240)
        par(mfcol=c(3,3))
        for (j in 1:ncol(dge$counts)) {
            xlim=c(0,quantile(dge$counts[,j],probs=.9,na.rm=T))
            plot(density(dge$counts[,j],na.rm=T),xlim=xlim,main=rownames(dge$samples)[j],col=colSam[j])
        }
        dev.off()

        png(paste("histogram_normCount",subsetName,".png",sep=""),width=3*240, height=3*240)
        par(mfcol=c(3,3))
        for (j in 1:ncol(dge$counts)) {
            xlim=c(0,quantile(dge$counts[,j],probs=.9,na.rm=T))
            xlim=c(0,quantile(dge$counts[,j],probs=.75,na.rm=T))
            hist(dge$counts[,j],xlim=xlim,main=rownames(dge$samples)[j])
        }
        dev.off()
    }

    dge$samples$group <- samples$treat
    
    grpUniq=unique(dge$samples$group)
    grpName=grpUniq
    
    for (grpId1 in 1:(length(grpUniq)-1)) {
        for (grpId2 in (grpId1+1):length(grpUniq)) {
            compFlag=paste("_",grpName[grpId2],"V",grpName[grpId1],sep="")
            
#			dgeT <- dge[, !is.na(dge$samples$group)]
            dgeT <- dge[,which(dge$samples$group%in%grpUniq)]
            id=rownames(dge1$samples)[which(dge1$samples$group%in%grpUniq[c(grpId1,grpId2)])]
            j=which(id%in%rownames(dgeT$samples))
            subsetName2=""
            subsetName2=subsetName
            
            dgeT <- estimateDisp(dgeT, trend = "none", robust = TRUE)
            print("table(is.infinite(dgeT$prior.df))")
            print(table(is.infinite(dgeT$prior.df)))
            summary(dgeT$tagwise.dispersion[is.infinite(dgeT$prior.df)])
            summary(dgeT$tagwise.dispersion[!is.infinite(dgeT$prior.df)])
            summary(dgeT$counts[is.infinite(dgeT$prior.df),])
            summary(dgeT$counts[!is.infinite(dgeT$prior.df),])
            summary(dgeT$counts[is.infinite(dgeT$prior.df),])

            "
            Dear Koen, 
            The results you are seeing are as one would expect for a technical datasets.
            Remember that the dispersion measures biological variation. The replicates
            here are just re-sequencing (and perhaps some re-preparation) of exactly the 
            same RNA samples, so there is no biological variation apart from slight 
            inaccuracies in the preparation of the samples. Hence the tagwise dispersions 
            will be very small (or even zero) and nearly equal. Hence the prior df will 
            be estimated to be large or even infinity. 
            Best wishes Gordon
            "

            summary(dgeT$prior.df)
            i=which(!is.infinite(dgeT$prior.df))
            summary(dgeT$prior.df[i])
            sqrt(dgeT$common.disp)

#			fName2=paste("_",organism,compFlag,subsetName2,sep="")
            fName2=paste(compFlag,subsetName2,sep="")
            save(dgeT,file=paste("dge_",organism,fName2,".RData",sep=""))

            ## -----------------------------------------


            if (organism%in%c("Homo_sapiens")) {
#				et <- exactTest(dgeT, pair = c("Wild type", "TGFbeta"))
                
#				dispersion - 
#				"common", "trended", "tagwise"
#				"auto" - Default behavior. Is to use most complex dispersions found in data object
                
                for (dispFlag in c("auto","common")) {
                    fName3=paste(fName2,"_",dispFlag,"Disp",sep="")
                    et <- exactTest(dgeT,dispersion=dispFlag,pair=grpUniq[c(grpId1,grpId2)])

                    if (F) {
                        source("funcs.R")

                        top <- topTags(et, n = 100)
                        top <- biomartify(top)

                        top <- topTags(et, n = nrow(et))
                        top <- biomartify(top, organism = organism)
                        write.table(top,paste("stat2",fName3,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)
#							save.image(paste("tmp",fName3,".RData",sep=""))
                    }

                    source("biomartify.1.2.R")
                    top <- topTags(et, n = 100)
                    top <- biomartify(top, organism = ifelse(organism=="Mus_musculus","MusMusculus","HomoSapiens"))
                    
                    source("biomartify.1.2.R")
                    top <- topTags(et, n = nrow(et))
                    top <- biomartify(top, organism = ifelse(organism=="Mus_musculus","MusMusculus","HomoSapiens"))
                    write.table(top,paste("stat",fName3,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)
                }
                
                ## Voom
                
                x=apply(dgeT$counts,1,max,na.rm=T)
                png("tmp.png")
                hist(x)
                dev.off()
                
                #for (minCnt in c(2,5,10)) {
                for (minCnt in c(0,20)) {
                    dgeF=dgeT[which(x>=minCnt),]
                    fName3=paste(fName2,"_voom_minCnt",minCnt,sep="")
                    design=model.matrix(~group,data=dgeF$samples)
                    if (F) {
                        fit <- voom(dgeF,design,plot=TRUE)
                        fit <- voom(dgeF,design,save.plot=TRUE)
                        png("tmp.png")
                        plot(fit)
                        dev.off()
                    }
                    dat <- voom(dgeF,design,save.plot=F)
                    save(dat,design,file=paste("voom",fName3,".RData",sep=""))
                    fit <- lmFit(dat,design)
                    rm(dat)
                    fit <- eBayes(fit)
                    #topTable(fit,coef=ncol(design))
                    colId=2
                    top=cbind(geneId=rownames(fit$coef),logFC=fit$coef[,colId],PValue=fit$p.value[,colId])
                    write.table(top,paste("stat",fName3,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)
                }
            }
            
            ## -----------------------------------------
            if (subsetFlag=="") {
                ## NOT USED
                source("biomartify.1.2.R")
                y=apply(dgeT$counts,1,function(x) {
                    mean(log2(x),na.rm=T)
                })
                ann <- biomartify(data.frame(geneId=rownames(dgeT$counts),logFC=y,logCPM=y,PValue=y,FDR=y,stringsAsFactors=F), organism = ifelse(organism=="Mus_musculus","MusMusculus","HomoSapiens"))
                ann <- ann[,which(!names(ann)%in%c("logFC","logCPM","PValue","FDR"))]
                write.table(ann,paste("ann_",organism,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)
            }
        }
    }
    
    ## -----------------------------------------
    ## Add logCPM for each group to output
    
    dispFlag <- "auto"
    #fName3 <- paste(fName2, "_", dispFlag, "Disp", sep = "")
    switch(cohort,
        "mmu"={
            varId="geneId"
            fName3="_transVsWt"
        },
        "mmu_hg"={
            varId="geneSymbol"
            fName3="_humanVsMouse"
        }
    )
    top <- read.table(paste("results/comparison/tmp/stat",fName3,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
    grpName2=rep("",length(grpUniq))
    grpName2[which(grpUniq=="Wild type")]="_wt"
    grpName2[which(grpUniq=="TGFbeta")]="_trans"
    out=matrix(nrow=nrow(top),ncol=length(grpUniq),dimnames=list(rownames(dge),paste("logCPM",grpName2,sep="")))
    for (grpId1 in 1:length(grpUniq)) {
        dgeT <- dge[,which(dge$samples$group%in%grpUniq[grpId1])]
        dgeT <- estimateDisp(dgeT, trend = "none", robust = TRUE)
        out[,grpId1]=dgeT$AveLogCPM
    }
    k=grep("foldChange",names(top))[1]
    top=cbind(top[,1:(k-1)],out[match(toupper(top[,varId]),toupper(rownames(out))),],top[,k:ncol(top)])
    write.table(top,paste("stat",fName3,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)
}

## -----------------------------------------
