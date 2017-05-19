library("Biobase")
library("viper")
args <- commandArgs(trailingOnly=TRUE)

logFC <- read.csv(args[1], row.names='Entrez_Id', stringsAsFactors=FALSE)

x <- logFC[,1]
names(x) <- rownames(logFC)

# brca_regulon <- load('viper/regulons/brca-tf-regulon.rda')
regulon <- load(args[2])

mrs <- msviper(x, regulon[1], verbose=FALSE)
result <- summary(mrs, mrs=length(mrs$es$nes))
write.csv(result, args[3])