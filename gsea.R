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
    library(biomaRt)
    listMarts()

    do:mart = useMart("ensembl")
    listDatasets(mart)
}



library(limma)
library(qvalue)


datadir="results/final/misc/"
datadir=""
ann <- read.table(paste(datadir,"ann_Homo_sapiens.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)

datadir="docs/gsea/"
refDBList=dir(datadir,pattern="rdata")
for (k in 1:length(refDBList)) {
    load(paste(datadir,refDBList[k],sep=""))
}


grpUniq="SURVIVAL"
grpUniq="APOPTOSIS"
grpUniq="TGF"
refDBName=rep("",length(refDBList))
for (dsRef in c("Hs.H","Hs.c1","Hs.c2","Hs.c3","Hs.c4","Hs.c6","Hs.c7")) {
    cat("\n\n##########################################\nRef geneset database: ",dsRef,"\n##########################################\n",sep="")
    switch(dsRef,
    "Hs.H"={
        refDBName[which(sub("human_","Hs.",sub("_v5p1.rdata","",refDBList))==dsRef)]="H hallmark gene sets"
        idGs <- ids2indices(Hs.H, ann$entrezgene[iA])
    },
    "Hs.c1"={
        refDBName[which(sub("human_","Hs.",sub("_v5p1.rdata","",refDBList))==dsRef)]="C1 positional gene sets"
        idGs <- ids2indices(Hs.c1, ann$entrezgene[iA])
    },
    "Hs.c2"={
        refDBName[which(sub("human_","Hs.",sub("_v5p1.rdata","",refDBList))==dsRef)]="C2 curated gene sets"
        idGs <- ids2indices(Hs.c2, ann$entrezgene[iA])
    },
    "Hs.c3"={
        refDBName[which(sub("human_","Hs.",sub("_v5p1.rdata","",refDBList))==dsRef)]="C3 motif gene sets"
        idGs <- ids2indices(Hs.c3, ann$entrezgene[iA])
    },
    "Hs.c4"={
        refDBName[which(sub("human_","Hs.",sub("_v5p1.rdata","",refDBList))==dsRef)]="C4 computational gene sets"
        idGs <- ids2indices(Hs.c4, ann$entrezgene[iA])
    },
    "Hs.c6"={
        refDBName[which(sub("human_","Hs.",sub("_v5p1.rdata","",refDBList))==dsRef)]="C6 oncogenic signatures"
        idGs <- ids2indices(Hs.c6, ann$entrezgene[iA])
    },
    "Hs.c7"={
        refDBName[which(sub("human_","Hs.",sub("_v5p1.rdata","",refDBList))==dsRef)]="C7 immunologic signatures"
        idGs <- ids2indices(Hs.c7, ann$entrezgene[iA])
    }
    )
    #print(length(grep(grpUniq,names(idGs))))
    print(names(idGs)[1:5])
}

timePt=36
minCnt=10
for (dsRef in c("Hs.H","Hs.c1","Hs.c2","Hs.c3","Hs.c4","Hs.c6","Hs.c7")) {
    cat("\n\n##########################################\nRef geneset database: ",refDBName[which(sub("human_","Hs.",sub("_v5p1.rdata","",refDBList))==dsRef)],"\n##########################################\n",sep="")
    for (timePt in c(4,8,12,36)) {
        cat("\n\n==========================================\nTime point: ",timePt," hr\n==========================================\n",sep="")
        
        datadir=""
        datadir="results/final/tables/"
        stat <- read.table(paste(datadir,"stat_TGFbetaVuntreated_",timePt,"hrs.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
        iA=match(stat$geneId,ann$geneId)

        datadir="/Users/royr/Downloads/tmp/"
        datadir="/voom_tmp/"
        load(paste(datadir,"voom_TGFbetaVuntreated_",timePt,"hrs_voom_minCnt",minCnt,".RData",sep=""))
        dim(dat)
        dim(stat)
        table(stat$geneId%in%rownames(dat))
        table(rownames(dat)%in%stat$geneId)
        iV=match(stat$geneId,rownames(dat))
        if (any(is.na(iV))) cat("Missing genes from voom !!!")

        switch(dsRef,
        "Hs.H"={
            idGs <- ids2indices(Hs.H, ann$entrezgene[iA])
        },
        "Hs.c1"={
            idGs <- ids2indices(Hs.c1, ann$entrezgene[iA])
        },
        "Hs.c2"={
            idGs <- ids2indices(Hs.c2, ann$entrezgene[iA])
        },
        "Hs.c3"={
            idGs <- ids2indices(Hs.c3, ann$entrezgene[iA])
        },
        "Hs.c4"={
            idGs <- ids2indices(Hs.c4, ann$entrezgene[iA])
        },
        "Hs.c6"={
            idGs <- ids2indices(Hs.c6, ann$entrezgene[iA])
        },
        "Hs.c7"={
            idGs <- ids2indices(Hs.c7, ann$entrezgene[iA])
        }
        )

        nRot=99
        nRot=9999
        rPV=romer(y=dat,index=idGs,design=design,contrast=2,nrot=nRot)
        rAnn=data.frame(id=rownames(rPV),numGene=rPV[,"NGenes"],stringsAsFactors=F)
        rPV=rPV[,-1]
        rownames(rAnn)=rownames(rPV)=NULL
        colnames(rPV)=tolower(colnames(rPV))
        rQV=matrix(nrow=nrow(rPV),ncol=ncol(rPV))
        colnames(rQV)=colnames(rPV)
        #for (k in 1:ncol(rPV)) {
        for (k in 3) {
            i=which(!is.na(rPV[,k]))
            if (F) {
            res=try(qvalue(rPV[i,k]))
            if (class(res)=="try-error") {
                #rQV[i,k]=p.adjust(rPV[i,k],method="BH")
                res=qvalue(rPV[i,k],lambda = seq(0,0.05,.01))
                rQV[i,k]=res$qvalues
            } else {
                rQV[i,k]=res$qvalues
            }
            }
            rQV[i,k]=p.adjust(rPV[i,k],method="BH")
        }
        fName=paste("gsea_TGFbetaVuntreated_",timePt,"hrs_refMSigDb",dsRef,".RData",sep="")
        save(rAnn,rPV,rQV,file=fName)


        k=which(colnames(rPV)=="mixed")
        
        pThres=0.05
        cat("No. of genesets with q-value < ",pThres,": ",sum(rQV[,k]<pThres),"\n",sep="")
        i=order(rPV[,k])
        i=i[rQV[i,k]<pThres]
        if (length(i)!=0) {
            out=unlist(lapply(idGs[i],function(x) {
                k=match(x,ann$entrezgene[iA])
                y=unique(ann$hgnc_symbol[iA][k[!is.na(k)]])
                c(length(y),paste(y[y!=""],collapse=","))
            }))
            fName=paste("gsea_TGFbetaVuntreated_",timePt,"hrs_refMSigDb",dsRef,"_fdr",pThres,".txt",sep="")
            tbl=cbind(geneset=rAnn$id[i],numGeneInData=out[seq(1,length(out),by=2)],pv=rPV[i,k],fdr=round(rQV[i,k],2),genes=out[seq(2,length(out),by=2)])
            write.table(paste("Reference dataset from Molecular Signatures Database (MSigDB): ",refDBName[which(sub("human_","Hs.",sub("_v5p1.rdata","",refDBList))==dsRef)],sep=""), file=fName,col.names=F,row.names=F, sep="\t",quote=F)
            write.table(paste(timePt," hr: TGFbeta vs. untreated",sep=""), file=fName,col.names=F,row.names=F, sep="\t",quote=F,append=T)
            write.table(tbl, file=fName,col.names=T,row.names=F, sep="\t",quote=F,append=T)
        }
        #pThres=0.1
        #cat("No. of genesets with q-value < ",pThres,": ",sum(rQV[,k]<pThres),"\n",sep="")
        
        cat("Top 10 genesets:\n",sep="")
        #print(topRomer(r,alt="mixed"))
        i=order(rPV[,k])[1:10]
        print(cbind(rAnn[i,],pv=rPV[i,k],qv=rQV[i,k]))
    }
}

rAnn$id[grep("SURVIVAL",rAnn$id)]

x="AARS,ANGPT2,AR,ATRX,BAD,CD44,CHKB,ECI1,DOCK3,ENDOG,ESR2,ETS2,FLOT2,FN1,GGTA1P,HSPG2,IL1RAP,IL13RA1,FOXK2,INCENP,INPP4A,LY9,MAP1A,MARS,METTL1,MMP2,MMP17,MYO9A,NEO1,NFATC1,OMG,P4HB,PEX1,PPBP,PRKY,PSMD2,RAB27A,RBBP5,RBBP7,REL,RGS2,SHOX2,STIL,SLC9A3,SPARC,SPINK1,NEK4,TEF,TTF1,TYK2,YWHAH,LAGE3,HIST1H3A,DCHS1,EIF2S2,TIMELESS,MAP3K6,TGFBRAP1,CRIPT,ITM2A,IL27RA,ZNF518A,USP15,BCAP31,GDF11,BAIAP2,ANP32B,TACC2,RAI1,CD3EAP,AP3M2"
y=strsplit(x,",")[[1]]
y[grep("TGF",y)]


#####################################################################################
#####################################################################################
Question: Update to MSigDB derivates for romer()?
0
gravatar for Jon Manning
9 months ago by
Jon Manning • 20
United Kingdom

The files quoted in the romer() manual at http://bioinf.wehi.edu.au/software/MSigDB/ are very helpful, but out of step with the current v5 of MSigDB. Are there plans to update them, or should I crack on and make my own?

Many thanks,

Jon
Here is a Rscript I use to do that sort of thing. This is for human.
#####################################################################################


## usage:  R CMD Rscript parseBroad_human.R  thedir  wheretoput
## where thedir is the dir containing the broad entrez-gene based data (broad_by_eg)
## wheretoput is the dir that the output should be saved in - no need to create first.
## the args are space delimited

args <- commandArgs(TRUE)

if(length(args) != 2) stop(paste("\nUsage: R CMD Rscript parseBroad_human.R",
" <thedir> <wheretoput>\n",
"<thedir> = the directory containing Entrez Gene based Gene sets",
"(currently /data2/romer/broad_by_eg)\n",
"<wheretoput> = the directory in which to output the mapped gene sets\n\n",
"Note that these arguments are space-delimited, and the function is best run",
" under nohup. Example:\n",
"nohup R CMD Rscript parseBroad_human.R broad_by_eg human2 &\n\n", sep = ""),
call. = FALSE)

thedir <- args[1]
wheretoput <- args[2]


fls <- list.files(thedir, pattern = "gmt$", full.names = TRUE)
thenames <- c("Positional", "Chemical_genetic_perturbations","Biocarta",
"KEGG", "Reactome","Canonical_pathway","microRNA",
"Transcription_factor_targets","Cancer_gene_neighborhoods","Cancer_modules",
"GO_biological_process","GO_cellular_component", "GO_molecular_function",
"Oncogenic", "Immunologic", "Hallmark")

parseBroad <- function(thefile){
    len <- as.numeric(strsplit(system(paste("wc -l", thefile), intern=TRUE), " ")[[1]][1])
    con <- file(thefile, "r")
    
    lst <- list(length = len)
    nam <- vector(length = len)
    
    for(i in seq_len(len)){
        #browser()
        theline <- readLines(con, 1L)
        theline <- strsplit(theline, "\t")
        nam[i] <- theline[[1]][1]
        lst[[i]] <- theline[[1]][-c(1:2)]
    }
    close(con)
    names(lst) <- nam
    return(lst)
    
}

broad <- lapply(fls, parseBroad)
## the canonical pathways include reactome, biocarta and KEGG, so subset
ind <- !names(broad[[6]]) %in% c(names(broad[[3]]), names(broad[[4]]), names(broad[[5]]))
broad[[6]] <- broad[[6]][ind]

names(broad) <- thenames
save(broad, file = paste0(wheretoput, "/broad_human_byset.Rdata"))

#####################################################################################
#####################################################################################
if(require("org.Hs.eg.db")) alias2Symbol(c("PUMA","NOXA","BIM"))



# org.Hs.egPATH - Mappings between Entrez Gene identifiers and KEGG pathway identifiers
## select() interface:
## Objects in this package can be accessed using the select() interface
## from the AnnotationDbi package. See ?select for details.
## Bimap interface:
x <- org.Hs.egPATH
# Get the entrez gene identifiers that are mapped to a KEGG pathway ID
mapped_genes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_genes])
if(length(xx) > 0) {
    # Get the PATH for the first five genes
    xx[1:5]
    # Get the first one
    xx[[1]]
}
# For the reverse map:
# Convert the object to a list
xx <- as.list(org.Hs.egPATH2EG)
# Remove pathway identifiers that do not map to any entrez gene id
xx <- xx[!is.na(xx)]
if(length(xx) > 0){
    # The entrez gene identifiers for the first two elements of XX
    xx[1:2]
    # Get the first one
    xx[[1]]
}
#####################################################################################
#####################################################################################

k=3
for (dsRef in c("Hs.H","Hs.c1","Hs.c2","Hs.c3","Hs.c4","Hs.c6","Hs.c7")) {
    cat("\n\n##########################################\nRef geneset database: ",refDBName[which(sub("human_","Hs.",sub("_v5p1.rdata","",refDBList))==dsRef)],"\n##########################################\n",sep="")
    for (timePt in c(4,8,12,36)) {
        cat("\n\n==========================================\nTime point: ",timePt," hr\n==========================================\n",sep="")
        datadir="results/gsea/"
        datadir=""
        fName1=paste("gsea_TGFbetaVuntreated_",timePt,"hrs_refMSigDb",dsRef,sep="")
        load(file=paste(datadir,fName1,".RData",sep=""))
        if (F) {
            x=p.adjust(rPV[,k],method="BH")
            print(table(rPV[,k]<0.05 & x>=0.05))
            print(table(rPV[,k]>=0.05 & x<0.05))
            print(table(rPV[,k]>=0.05 & rQV[,k]<0.05))
            nm="Q-value"
        }
        rQV[,k]=p.adjust(rPV[,k],method="BH")
        nm="Adjusted p-value"
        png(paste(fName1,".png",sep=""))
        lim=c(0,1); lim=NULL
        plot(rPV[,k],rQV[,k],xlim=lim,ylim=lim,main=paste(dsRef,", ",timePt," hr: TGFbeta vs. untreated",sep=""),xlab="P-value",ylab=nm)
        abline(c(0,1))
        dev.off()
        

        pThres=0.05
        cat("No. of genesets with adjusted p-value < ",pThres,": ",sum(rQV[,k]<pThres),"\n",sep="")
        i=order(rPV[,k])
        i=i[rQV[i,k]<pThres]
        if (length(i)!=0) {
            out=unlist(lapply(idGs[i],function(x) {
                k=match(x,ann$entrezgene[iA])
                y=unique(ann$hgnc_symbol[iA][k[!is.na(k)]])
                c(length(y),paste(y[y!=""],collapse=","))
            }))
            fName=paste(fName1,"_fdr",pThres,".txt",sep="")
            tbl=cbind(geneset=rAnn$id[i],numGeneInData=out[seq(1,length(out),by=2)],pv=rPV[i,k],fdr=signif(rQV[i,k],2),genes=out[seq(2,length(out),by=2)])
            write.table(paste("Reference dataset from Molecular Signatures Database (MSigDB): ",refDBName[which(sub("human_","Hs.",sub("_v5p1.rdata","",refDBList))==dsRef)],sep=""), file=fName,col.names=F,row.names=F, sep="\t",quote=F)
            write.table(paste(timePt," hr: TGFbeta vs. untreated",sep=""), file=fName,col.names=F,row.names=F, sep="\t",quote=F,append=T)
            write.table(tbl, file=fName,col.names=T,row.names=F, sep="\t",quote=F,append=T)
        }
    }
}

#####################################################################################
#####################################################################################
