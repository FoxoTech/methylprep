# BiocManager::install(c("minfi","minfiData"))
library(minfi)
# specify directory
#minfi_baseDir = paste0("/Users/nrigby/Desktop/idats_standard/batch_1052641/")
baseDir = paste0("/Users/mmaxmeister/methylprep/docs/example_data/GSE69852/minfi/")
# read samplesheet
targets = read.metharray.sheet(baseDir)
# read IDAT's into RGChannelSet
rgSet <- read.metharray.exp(targets = targets, verbose = TRUE)
# preprocessRaw
mSet.raw = preprocessRaw(rgSet)
# make files
meth.raw = getMeth(mSet.raw)
write.csv(file=file.path(baseDir,'minfi_raw_meth.csv'),x=meth.raw,row.names=TRUE)
unmeth.raw = getUnmeth(mSet.raw)
write.csv(file=file.path(baseDir,'minfi_raw_unmeth.csv'),x=unmeth.raw,row.names=TRUE)
betas.raw = getBeta(mSet.raw)
write.csv(file=file.path(baseDir,'minfi_raw_betas.csv'),x=betas.raw,row.names=TRUE)
# preprocessNoob
mSet.noob = preprocessNoob(rgSet)
# make files
meth.noob = getMeth(mSet.noob)
write.csv(file='~/Desktop/idats_standard/minfi_noob_meth.csv',x=meth.noob,row.names=TRUE,col.names=TRUE)
unmeth.noob = getUnmeth(mSet.noob)
write.csv(file='~/Desktop/idats_standard/minfi_noob_unmeth.csv',x=unmeth.noob,row.names=TRUE,col.names=TRUE)
betas.noob = getBeta(mSet.noob)
write.csv(file='~/Desktop/idats_standard/minfi_noob_betas.csv',x=betas.noob,row.names=TRUE,col.names=TRUE)
