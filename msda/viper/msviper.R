library("Biobase")
library("viper")
args <- commandArgs(trailingOnly=TRUE)

exprsFile <- args[1]
exprs <- as.matrix(read.table(exprsFile, header=TRUE, sep=',', row.names='Gene_Symbol', as.is=TRUE))


pDataFile <- args[2]
pData <- read.table(pDataFile, row.names='Sample', header=TRUE, sep=',')

all(rownames(pData)==colnames(exprs))

phenoData <- new("AnnotatedDataFrame", data=pData)

exSet <- ExpressionSet(assayData=exprs, phenoData=phenoData)

exSet

signature <- rowTtest(exSet, args[3] , args[4], args[5])

# signature <- rowTtest(exSet, "Molecular_subtype", c("Basal", "Basal A", "Basal B"), "Luminal")

signature <- (qnorm(signature$p.value/2, lower.tail=FALSE) * sign(signature$statistic))[, 1]

# nullmodel <- ttestNull(exSet, "Molecular_subtype", c("Basal", "Basal A", "Basal B"), "Luminal", per=1000, repos=TRUE, verbose=FALSE)

nullmodel <- ttestNull(exSet, args[3], args[4], args[5], per=1000, repos=TRUE, verbose=FALSE)


chris_regulon <- load('viper/regulon_symbol.rdata')

mrs <- msviper(signature, regulon1, nullmodel, verbose=FALSE)
result <- summary(mrs, mrs=length(mrs$es$nes))
write.csv(result, args[6])