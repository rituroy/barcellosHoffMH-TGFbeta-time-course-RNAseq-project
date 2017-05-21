## Get DE genes for pairs of sample groups for RNAseq count data

R.utils::use("aroma.seq, edgeR, aroma.light, matrixStats")

verbose=T
verbose=F

dataset <- "tmp"
cohort <- "hg"
organism <- "Homo_sapiens"

counts <- HTSeqCountDataSet$byName(dataset, tags = "tophat2,gtf", organism = organism)
counts <- setFullNamesTranslator(counts, function(name, ...) chartr(";", "/", name))
counts

data <- extractMatrix(counts, column = 2L, colClass = "integer")
dge <- extractDGEList(counts)

db <- TabularTextFile("data/sampleInfo.txt")
samples <- readDataFrame(db)
samples <- samples[match(colnames(data), samples$id), ]
samples1 <- samples

## ------------
## Check correlation among samples

rho <- cor(data, method = "spearman")
round(summary(as.vector(rho)),2)
j=which(samples$time==4)
round(rho[j,j],2)

## ------------

dge$samples[, -1]
samples <- cbind(samples, dge$samples[, c("lib.size", "norm.factors")])
if (verbose) {
    write.table(cbind(geneId = rownames(dge$counts), dge$counts), file = paste("count_raw_", organism, ".txt", sep = ""), col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
    write.table(samples, file = paste("sample_", organism, ".txt", sep = ""), col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
}

dge <- calcNormFactors(dge, method = "TMM")

samples <- samples1[match(rownames(dge$samples), samples1$id), ]
samples <- cbind(samples, dge$samples[, c("lib.size", "norm.factors")])
if (verbose) {
    write.table(samples, file = paste("sample_", organism, ".txt", sep = ""), col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
}
dge$samples$group <- samples$treat

dgeAll <- dge


#--------------------------
dge=dgeAll

samplesThis=samples[match(colnames(dge),samples$id),]
samplesThis$treat2=2-as.integer(as.factor(samplesThis$treat))
j=1:ncol(dge)
j=j[which(!samplesThis$id[j]%in%c("Donor_1_TGFbeta_4hrs","Donor_3_4_hrs"))]
dge=dge[,j]
samplesThis=samplesThis[j,]
j=order(samplesThis$donor,samplesThis$treat2,samplesThis$time)
dge=dge[,j]
samplesThis=samplesThis[j,]
lcpmT=cpm(dge,log=TRUE)
tbl=as.data.frame(lcpmT)
tbl=cbind(ensembl_gene_id=rownames(tbl),tbl)
if (verbose) {
    write.table(tbl, file=paste("logCPM.txt",sep=""),col.names=T,row.names=F, sep="\t",quote=F)
    gzip("logCPM.txt")
}
samplesThis=samples[match(colnames(dge),samples$id),]
samplesThis$treat2=2-as.integer(as.factor(samplesThis$treat))
j=1:ncol(dge)
j=j[which(!samplesThis$id[j]%in%c("Donor_1_TGFbeta_4hrs","Donor_3_4_hrs"))]
dge=dge[,j]
samplesThis=samplesThis[j,]
j=order(samplesThis$donor,samplesThis$treat2,samplesThis$time)
dge=dge[,j]
samplesThis=samplesThis[j,]
tbl=as.data.frame(dge$counts)
tbl=cbind(ensembl_gene_id=rownames(tbl),tbl)
if (verbose) {
    write.table(tbl, file=paste("countRaw.txt",sep=""),col.names=T,row.names=F, sep="\t",quote=F)
    gzip("countRaw.txt")
}
#--------------------------

dge=dgeAll


## -----------------------------------------

subsetFlag <- "_4hrs"
subsetFlag=""

subsetList=paste("_",sort(unique(samples1$time)),"hrs",sep="")

for (subsetFlag in subsetList) {

    fName1 <- paste("_", organism, subsetFlag, sep = "")
    dge <- dgeAll
    samples <- samples1[match(rownames(dge$samples), samples1$id), ]
    samples <- cbind(samples, dge$samples[, c("lib.size", "norm.factors")])

    j <- 1:nrow(dge$samples)
    if (subsetFlag!="") j=which(paste("_",samples$time,"hrs",sep="")==subsetFlag)
    j=j[which(!rownames(dge$samples)[j]%in%c("Donor_1_TGFbeta_4hrs","Donor_3_4_hrs"))]
    subsetName=subsetFlag
    dge <- dge[, j]
    samples <- samples[j,]
    print(subsetFlag)
    print(rownames(dge$samples))

    colSam <- rainbow(ncol(dge$counts))


    dge <- calcNormFactors(dge, method = "TMM")

    dge$samples$group <- samples$treat

    grpUniq <- unique(dge$samples$group)
    grpName <- grpUniq

    ## Loop to perform association test between pairs of sample groups
    for (grpId1 in 1:(length(grpUniq)-1)) {
        for (grpId2 in (grpId1+1):length(grpUniq)) {
            compFlag <- paste("_", grpName[grpId2], "V", grpName[grpId1], sep = "")
            
            samId=which(dge$samples$group%in%grpUniq)
            dgeT <- dge[,samId]
            id <- rownames(dgeAll$samples)[which(dgeAll$samples$group%in%grpUniq[c(grpId1, grpId2)])]
            j <- which(id%in%rownames(dgeT$samples))
            subsetName2 <- subsetName
            
            dgeT <- estimateDisp(dgeT, trend = "none", robust = TRUE)
            print("table(is.infinite(dgeT$prior.df))")
            print(table(is.infinite(dgeT$prior.df)))
            summary(dgeT$tagwise.dispersion[is.infinite(dgeT$prior.df)])
            summary(dgeT$tagwise.dispersion[!is.infinite(dgeT$prior.df)])
            summary(dgeT$counts[is.infinite(dgeT$prior.df), ])
            summary(dgeT$counts[!is.infinite(dgeT$prior.df), ])
            summary(dgeT$counts[is.infinite(dgeT$prior.df), ])

            summary(dgeT$prior.df)
            i <- which(!is.infinite(dgeT$prior.df))
            summary(dgeT$prior.df[i])
            sqrt(dgeT$common.disp)

            fName2 <- paste(compFlag, subsetName2, sep = "")
            if (verbose) {
                save(dgeT, file = paste("dge_", organism, fName2, ".RData", sep = ""))
            }

            #cntVec=apply(dgeT$counts,1,max,na.rm=T)
            cntVec=apply(t(t(dgeT$counts)*dgeT$sample$norm.factors),1,max,na.rm=T)
            
            ## -----------------------------------------
            if (F) {
                dispFlag <- "common"
                dispFlag <- "tagwise"
                dispFlag <- "trended"
                dispFlag <- "auto"
                for (minCnt in c(0,2,5,10,20)) {
                    fName3 <- paste(fName2, "_", dispFlag, "Disp_minCnt",minCnt,sep = "")
                    dgeF=dgeT[which(cntVec>=minCnt),]
                    dgeF <- estimateDisp(dgeF, trend = "none", robust = TRUE)
                    et <- exactTest(dgeF, dispersion = dispFlag, pair = grpUniq[c(grpId1, grpId2)])
                    
                    #source(paste(dirSrc,"functions/biomartify.1.1.R",sep=""))
                    source(paste(dirSrc,"functions/biomartify.1.2.R",sep=""))
                    top <- topTags(et, n = nrow(et))
                    write.table(cbind(geneId=rownames(top),top$table), paste("stat", fName3, ".txt", sep = ""), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
                    #top <- biomartify(top, organism = "HomoSapiens")
                    #write.table(top, paste("stat", fName3, ".txt", sep = ""), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
                    
                    if (F) {
                        source(paste(dirSrc,"functions/biomartify.1.2.R",sep=""))
                        i=1:10
                        top1 <- biomartify(top[i,], organism = "HomoSapiens")
                        top1 <- biomartify(as.data.frame(top[i,]), organism = "HomoSapiens")
                        
                        dataset <- "hsapiens_gene_ensembl"
                        ensembl <- useMart("ensembl", dataset = dataset)
                        attrs <- c("ensembl_gene_id", "hgnc_symbol", "hgnc_id", "chromosome_name", "start_position", "end_position",
                        "strand", "gene_biotype", "description")
                        res=getBM(attributes=attrs, filters = "ensembl_gene_id", values = top$geneId[i], mart=ensembl, curl = NULL, checkFilters = TRUE, verbose = FALSE, uniqueRows = TRUE, bmHeader = FALSE)

                        top1 <- biomartify(top, organism = "HomoSapiens")
                        write.table(top1, paste("stat", fName3, ".txt", sep = ""), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
                    }
                }
            }
            
            ## -----------------------------------------
            if (T) {
                effectFlag="random"
                effectFlag="fixed"
                if (effectFlag=="random" & subsetFlag=="_4hrs") break
                minCntList=c(0,2,5,10,20)
                minCntList=c(10)
                for (minCnt in minCntList) {
                    fName3=paste(fName2,"_voom_minCnt",minCnt,sep="")
                    dgeF=dgeT[which(cntVec>=minCnt),]
                    dgeF <- estimateDisp(dgeF, trend = "none", robust = TRUE)
                    group=2-as.integer(as.factor(dgeF$samples$group))
                    if (F) {
                        fit <- voom(dgeF,design,plot=TRUE)
                        fit <- voom(dgeF,design,save.plot=TRUE)
                        png("tmp.png")
                        plot(fit)
                        dev.off()
                    }
                    if (effectFlag=="fixed") {
                        design=model.matrix(~group+samples$donor[samId])
                        objThis <- voom(dgeF,design,save.plot=F)
                        fit <- lmFit(objThis,design)
                    } else {
                        design=model.matrix(~group)
                        objThis <- voom(dgeF,design,save.plot=F)
                        corfit <- duplicateCorrelation(objThis,design,block=samples$donor[samId])
                        cat("corfit$consensus: ",corfit$consensus,sep="")
                        fit <- lmFit(objThis,design,block=samples$donor[samId],correlation=corfit$consensus)
                    }
                    dat <- objThis
                    save(dat,design,file=paste("voom",fName3,".RData",sep=""))
                    rm(dat)
                    fit <- eBayes(fit)
                    #topTable(fit,coef=ncol(design))
                    colId=2
                    #top=data.frame(geneId=rownames(fit$coef),logFC=fit$coef[,colId],logCPM=dgeF$AveLogCPM,PValue=fit$p.value[,colId])
                    top=data.frame(geneId=rownames(fit$coef),logFC=fit$coef[,colId],logCPM=dgeF$AveLogCPM,t=fit$t[,colId],PValue=fit$p.value[,colId])
                    top$FDR=NA
                    library(qvalue)
                    i=which(!is.na(top$PValue))
                    top$FDR[i]=qvalue(top$PValue[i])$qvalues
                    write.table(top,paste("stat",fName3,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)
                    
                    if (F) {
                        top1 <- biomartify(top, organism = "HomoSapiens")
                        write.table(top1,paste("stat",fName3,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)
                    }
                }
            }

        }
    }
}

## -----------------------------------------
