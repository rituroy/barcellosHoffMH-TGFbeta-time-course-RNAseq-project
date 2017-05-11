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


##############################################
## Create sample info table

dirClin="docs/"
fNameClin="150903_D00372_0318_AC7FNUANXX_Combo_HSQ_98.txt"
phen=read.table(paste(dirClin,fNameClin,sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
names(phen)[match(c("Lane","Sample.ID","Sample.Ref","Index","Description","Control","Project","Yield..Mbases.","X..PF","X..Reads","X..of.raw.clusters.per.lane","X..Perfect.Index.Reads","X..One.Mismatch.Reads..Index.","X..of....Q30.Bases..PF.","Mean.Quality.Score..PF."),names(phen))]=c("lane","id","sampleRef","index","desc","control","project","yield","percPF","numReads","percRawClustersPerLane","percPerfectIndexReads","percOneMismatchReads","percQ30Bases","meanQualityScore")

phen=phen[grep("Donor",phen$id),]

phen$numReads=gsub(",","",phen$numReads)
for (k in 1:ncol(phen)) {
    phen[,k]=sub(" +$", "", phen[,k])
}
for (k in which(names(phen)%in%c("lane","yield","percPF","numReads","percPerfectIndexReads","percOneMismatchReads"))) {
    phen[,k]=as.integer(phen[,k])
}
for (k in which(names(phen)%in%c("percRawClustersPerLane","percQ30Bases","meanQualityScore"))) {
    phen[,k]=as.numeric(phen[,k])
}

out=as.data.frame(t(sapply(phen$id,function(x) {
    y=gsub("Donor_|hrs","",sub("_hrs","hrs",x))
    y=strsplit(y,"_")[[1]]
    if (length(y)==2) y=c(y[1],"untreated",y[2])
    y
},USE.NAMES=F)),stringsAsFactors=F)
names(out)=c("donor","treat","time")
out$time=as.integer(out$time)
phen=cbind(phen,out)

write.table(phen, file="sampleInfo.txt", sep="\t", col.names=T, row.names=F, quote=F)

##############################################
## Check numbers from manuscript

colGeneId="geneId"; pThres=0.05
datadir="results/final/tables/"
dat1=read.table(paste(datadir,"stat_TGFbetaVuntreated_8hrs.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
#stat1=statVF_8_10
ann_2=ann10[match(dat1[,colGeneId],ann10[,colGeneId]),]
i2=which(dat1$Qvalue<pThres)
i11=which(!duplicated(dat1$hgnc_symbol) & dat1$hgnc_symbol!="")
i21=which(!duplicated(dat1$hgnc_symbol[i2]) & dat1$hgnc_symbol[i2]!="")
nrow(dat1) # Total no. of genes tested
length(i11) # Total no. of unique genes tested
length(i2) # No. of significant genes
length(i21) # No. of unique significant genes




##############################################
