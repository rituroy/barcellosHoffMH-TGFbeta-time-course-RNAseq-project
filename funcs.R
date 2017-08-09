getAnnotation=function(build=38) {
    phen=read.table(paste("docs/sampleInfo.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
    datadir="results/final/misc/"
    ann=read.table(paste(datadir,"ann_Homo_sapiens.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
    if (build==37) {
        library(biomaRt)

        dataset <- "hsapiens_gene_ensembl"
        #ensembl <- useMart("ensembl", dataset = dataset)
        ensembl = useEnsembl(biomart="ensembl", dataset=dataset, GRCh=37)
        #attrs <- c("ensembl_gene_id", "hgnc_symbol", "hgnc_id", "chromosome_name", "start_position", "end_position","strand", "gene_biotype", "description")
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
    }
    names(ann)[match(c("hgnc_symbol"),names(ann))]=c("geneSym")
    ann$chr=as.integer(ann$chromosome_name)
    ann$chr[which(ann$chromosome_name=="X")]=23
    ann$chr[which(ann$chromosome_name=="Y")]=24
    ann
}


#####################################################################################
#####################################################################################
