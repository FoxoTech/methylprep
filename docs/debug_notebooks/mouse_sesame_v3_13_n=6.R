library(BiocManager)
library(sesame)
library(preprocessCore)
# after upgrading R3.6-->R4.1 I had to move cache: https://bioconductor.org/packages/devel/bioc/vignettes/ExperimentHub/inst/doc/ExperimentHub.html#default-caching-location-update
# sesameDataCache("MM285") # run this once, before you can read mouse
# version 3.13 DOCS HERE: https://bioconductor.org/packages/devel/bioc/vignettes/sesame/inst/doc/sesame.html

in_dir = paste0('/Volumes/LEGX/Barnes/mouse_test/')

savefiles = function(sdfs, idx, thisname) {
  x_thisframe <- data.frame(rbind( IG(sdfs[[idx]]), IR(sdfs[[idx]]), II(sdfs[[idx]]) ))
  x_filename <- paste0('sesame_mouse_', thisname ,'.csv', sep="")
  x_file <- file.path(in_dir, x_filename)
  write.csv(x_file, x=x_thisframe, row.names=TRUE, na="")
}

sdfs <- lapply(searchIDATprefixes(in_dir), readIDATpair) #readIdat is for one file; readIDATpair is for pairs.
savefiles(sdfs, 1, 'raw')
# INFER
sdfs = lapply(sdfs, inferTypeIChannel, verbose=TRUE)
savefiles(sdfs, 1, 'infer')

# POOBAH ssets <- lapply(ssets, detectionMask) # poobah and qualityMask now set in sdf@extra$mask by readIDATpair
sdfs <- lapply(sdfs, pOOBAH)
savefiles(sdfs, 1, 'poobah')

### -- these don't return a modified sdf. so they need a new location.
# NOOB
sdfs.noob <- lapply(sdfs, noob) # creates a bunch more sdf objects in a list here
savefiles(sdfs.noob, 1, 'noob')

# DYE-BIAS
sdfs.dye <- lapply(sdfs.noob, dyeBiasCorrTypeINorm)
savefiles(sdfs.dye, 1, 'dye')

# BETAS
betas <- lapply(sdfs.dye, getBetas)
savefiles(betas, 1, 'beta')

Rprof(tf <- "log.log", memory.profiling = TRUE)
idat_dir = paste0('/Volumes/LEGX/Barnes/mouse_test/')
betas= openSesame(idat_dir, BPPARAM=BiocParallel::MulticoreParam(2))
Rprof (NULL) ; print(summaryRprof(tf))
system.time( openSesame(idat_dir, BPPARAM=BiocParallel::MulticoreParam(2)) )
# sesame_version user     system elapsed
# 3.13           42.691   2.705  23.291
