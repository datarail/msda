library("Biobase")
library("viper")
args <- commandArgs(trailingOnly=TRUE)

identifier <- args[2]
exprsFile <- args[1]
exprs <- as.matrix(read.table(exprsFile, header=TRUE, sep=',', row.names=identifier, as.is=TRUE))


regulon <- load(args[3])

vpres <- viper(exprs, get(regulon), verbose=FALSE)
write.csv(vpres, args[4])