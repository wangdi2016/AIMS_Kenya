### read in RNASeq table

##
## read gene expression and clean gene list
##
D <- read.table("./read_counts_TPM.txt", header=T, stringsAsFactors=F)
g <- read.table("./kenya.geneid.4AIMS.clean.tsv", header=T, stringsAsFactors=F)

## put genename into rowname
rownames(D) <- D[,1]

## copy gene_id column to genename so we can use merge later
D$genename <- D$gene_id

## merge two dataframes by genename  
df <- merge(D,g, by="genename")
write.table(df, file="kenya.count.4AIMS.tsv", sep="\t", quote=F)
dim(df)

## split df.new into Dx
Dx <- as.matrix(df[,3:47])

### read in EntrezID as vector
ID<-as.character(df$geneid)

### load AIMS
library(AIMS)
subtypes<-applyAIMS(Dx,ID)
head(subtypes$cl)
table(subtypes$cl)

write.table(file="kenya.subtypes.tpm.aims.txt", subtypes$cl, sep="\t", quote=F)
