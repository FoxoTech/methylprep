if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.13")

in_dir = paste0('/Volumes/LEGX/GEO/test_pipeline/GSE147430/GPL21145/sesame/')

sesameDataCache("EPIC")

savefiles = function(sdfs, idx, thisname) {
  x_thisframe <- data.frame(rbind( IG(sdfs[[idx]]), IR(sdfs[[idx]]), II(sdfs[[idx]]) ))
  x_filename <- paste0('sesame_mouse_', thisname ,'.csv', sep="")
  x_file <- file.path(in_dir, x_filename)
  write.csv(x_file, x=x_thisframe, row.names=TRUE, na="")
}

# RAW INTENSITY
sdfs <- lapply(searchIDATprefixes(in_dir), readIDATpair) #readIdat is for one file; readIDATpair is for pairs.
savefiles(sdfs, 1, 'raw')

# INFER
sdfs.infer = lapply(sdfs, inferTypeIChannel, verbose=TRUE)
savefiles(sdfs, 1, 'infer')

# POOBAH ssets <- lapply(ssets, detectionMask) # poobah and qualityMask now set in sdf@extra$mask by readIDATpair
sdfs <- lapply(sdfs, pOOBAH)
savefiles(sdfs, 1, 'poobah') # doesn't save the poobah filtered part

### -- these don't return a modified sdf. so they need a new location.
# NOOB
sdfs.noob <- lapply(sdfs, noob) # creates a bunch more sdf objects in a list here
savefiles(sdfs.noob, 1, 'noob') # noob doesn't actually exclude poobah probes here? unless auto

# DYE-BIAS
sdfs.dye <- lapply(sdfs.noob, dyeBiasCorrTypeINorm)
savefiles(sdfs.dye, 1, 'dye')

# BETAS
betas <- lapply(sdfs.dye, getBetas) # this uses @extra$mask
write.csv(file.path(in_dir, 'sesame_mouse_betas.csv'), x=betas, row.names=TRUE, na="")

