library("Biobase")
library("viper")
library(stringr)
args <- commandArgs(trailingOnly=TRUE)

exprsFile <- args[1]
identifier <- args[6]
exprs <- as.matrix(read.table(exprsFile, header=TRUE, sep=',', row.names=identifier, as.is=TRUE))

pDataFile <- args[2]
pData <- read.table(pDataFile, row.names='Sample', header=TRUE, sep=',')
all(rownames(pData)==colnames(exprs))
phenoData <- new("AnnotatedDataFrame", data=pData)

exSet <- ExpressionSet(assayData=exprs, phenoData=phenoData)

test_str <- args[4]
test_samples <- unlist(str_split(test_str, ","))
ref_str <- args[5]
ref_samples <- unlist(str_split(ref_str, ", "))
signature <- rowTtest(exSet, args[3], test_samples, ref_samples)

signature <- (qnorm(signature$p.value/2, lower.tail=FALSE) * sign(signature$statistic))[, 1]
# nullmodel <- ttestNull(exSet, args[3], test_samples, ref_samples, per=1000, repos=TRUE, verbose=FALSE)

regulon <- load(args[7])

mrs <- msviper(signature, get(regulon),  verbose=FALSE)
result <- summary(mrs, mrs=length(mrs$es$nes))
write.csv(result, args[8])