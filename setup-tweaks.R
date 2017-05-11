##############################################################
# Workound
# Currently aroma.seq assumes annotation data files to
# be in directory annotationData/organisms/<organism>/
# At CBC, we are keeping ours in subdirectories of this,
# which confuses aroma.seq.  For now, lets create links
# such that aroma.seq only sees the above structure.
##############################################################
library("R.utils")
library("aroma.seq")

pathS <- "/cbc/data/annotationData/organisms"
pathS <- Arguments$getReadablePath(pathS)

path <- "annotationData/organisms"
path <- Arguments$getWritablePath(path)


organism <- "Homo_sapiens"
assembly <- c(GRC="GRCh38", UCSC="hg38")
release <- 79

pathA <- file.path(path, organism)
target <- file.path(pathS, organism, fullname(assembly), "Ensembl", release)
createLink(pathA, target=target)

# Sanity check
filename <- sprintf("%s.%s.dna.chr=1-MT.fa", organism, assembly["GRC"])
fa <- FastaReferenceFile(filename, path=pathA)
print(fa)


organism <- "Mus_musculus"
assembly <- c(GRC="GRCm38", UCSC="mm10")
release <- 79

pathA <- file.path(path, organism)
target <- file.path(pathS, organism, fullname(assembly), "Ensembl", release)
createLink(pathA, target=target)

# Sanity check
filename <- sprintf("%s.%s.dna.chr=1-MT.fa", organism, assembly["GRC"])
fa <- FastaReferenceFile(filename, path=pathA)
print(fa)
