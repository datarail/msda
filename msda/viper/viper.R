library("Biobase")
library("viper")
args <- commandArgs(trailingOnly=TRUE)

chris_regulon <- load('viper/regulon_symbol.rdata')

exprsFile <- args[1]
exprs <- as.matrix(read.table(exprsFile, header=TRUE, sep=',', row.names='Gene_Symbol', as.is=TRUE))


pDataFile <- args[2]
pData <- read.table(pDataFile, row.names='Sample', header=TRUE, sep=',')

all(rownames(pData)==colnames(exprs))

phenoData <- new("AnnotatedDataFrame", data=pData)

exSet <- ExpressionSet(assayData=exprs, phenoData=phenoData)

exSet

vpres <- viper(exprs, regulon1, verbose=FALSE)

write.csv(vpres, args[3])