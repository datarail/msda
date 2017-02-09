library("Biobase")
library("viper")
args <- commandArgs(trailingOnly=TRUE)


logFC <- read.csv(args[1], row.names='Gene_Symbol', stringsAsFactors=FALSE)

x <- logFC[,1]
names(x) <- rownames(logFC)

brca_regulon <- load('viper/regulons/regul_symbol.rda')

mrs <- msviper(x, regul_symbol,  verbose=FALSE)
result <- summary(mrs, mrs=length(mrs$es$nes))
write.csv(result, args[2])