### read in RNASeq table

##
## read gene expression and clean gene list
##
D <- read.table("./read_counts_TPM.txt", header=T, stringsAsFactors=F)
g <- read.table("./kenya.geneid.4AIMS.clean.tsv", header=T, stringsAsFactors=F)

## split D into Dx number only
rownames(D) <- D[,1]

## copy gene_id column to genename so we can use merge later
D$genename <- D$gene_id

## merge two dataframes by genename  
df.new <- merge(D,g, by="genename")
dim(df.new)

Dx <- as.matrix(df.new[,3:47])

write.table(df.new, file="kenya.count.4AIMS.tsv", sep="\t", quote=F)

### read in EntrezID as vector
ID<-as.character(df.new$geneid)

### load AIMS
library(AIMS)
subtypes<-applyAIMS(Dx,ID)
head(subtypes$cl)
table(subtypes$cl)

write.table(file="kenya.subtypes.tpm.aims.txt", subtypes$cl, sep="\t", quote=F)
