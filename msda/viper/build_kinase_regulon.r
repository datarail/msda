args = commandArgs(trailingOnly=TRUE)

get_regulon <- function(regul_file){
  regul_mat <- read.delim(regul_file,sep = ',',as.is = T)
  regulon <- list()
  for(i in unique(regul_mat$KINASE)){
    regul_mat_i <- regul_mat[regul_mat$KINASE == i,]
    regulon[[i]] <- list('tfmode' = rep(1,nrow(regul_mat_i)),'likelihood'=regul_mat_i$confidence)
    names(regulon[[i]]$tfmode) <- paste(regul_mat_i$Gene_Symbol, regul_mat_i$Site, sep = '_')
  }
  return(regulon)
}


regulon_file <- args[1]
regulon <- get_regulon(regulon_file)
save(regulon, file = args[2])

## Command line example
## rscript filename.csv regulon.rdata


                                    
