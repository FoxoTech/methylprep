library(minfi)
# specify directory
minfi_baseDir = paste0("/Users/nrigby/Desktop/idats_standard/batch_1052641/")
# read samplesheet
minfi_targets = read.metharray.sheet(minfi_baseDir)
# read IDAT's into RGChannelSet
rgSet <- read.metharray.exp(targets = minfi_targets1,verbose = TRUE)
# preprocessRaw
mSet.raw = preprocessRaw(rgSet)
# make files
meth.raw = getMeth(mSet.raw)
write.csv(file='~/Desktop/idats_standard/minfi_raw_meth.csv',x=meth.raw1,row.names=TRUE,col.names=TRUE)
unmeth.raw = getUnmeth(mSet.raw)
write.csv(file='~/Desktop/idats_standard/minfi_raw_unmeth.csv',x=unmeth.raw1,row.names=TRUE,col.names=TRUE)
betas.raw = getBeta(mSet.raw)
write.csv(file='~/Desktop/idats_standard/minfi_raw_betas.csv',x=betas.raw1,row.names=TRUE,col.names=TRUE)
# preprocessNoob
mSet.noob = preprocessNoob(rgSet)
# make files
meth.noob = getMeth(mSet.noob)
write.csv(file='~/Desktop/idats_standard/minfi_noob_meth.csv',x=meth.noob,row.names=TRUE,col.names=TRUE)
unmeth.noob = getUnmeth(mSet.noob)
write.csv(file='~/Desktop/idats_standard/minfi_noob_unmeth.csv',x=unmeth.noob,row.names=TRUE,col.names=TRUE)
betas.noob = getBeta(mSet.noob)
write.csv(file='~/Desktop/idats_standard/minfi_noob_betas.csv',x=betas.noob,row.names=TRUE,col.names=TRUE)
