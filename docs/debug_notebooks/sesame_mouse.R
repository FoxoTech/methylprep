library(BiocManager)
library(sesame)
library(preprocessCore)
in_dir = paste0('/Users/mmaxmeister/methylprep/docs/example_data/mouse/')

ssets <- lapply(searchIDATprefixes(in_dir), readIDATpair)
ssets.noob <- lapply(ssets, noob)
ssets <- lapply(ssets, detectionMask) # poobah
ssets.dye <- lapply(ssets.noob, dyeBiasCorrTypeINorm)
ssets.noob.df = data.frame(rbind( IG(ssets.noob[[1]]), IR(ssets.noob[[1]]), II(ssets.noob[[1]]) ))
write.csv(file= file.path(in_dir, 'sesame_mouse_noob.csv'), x=ssets.noob.df, row.names=TRUE)
ssets.dye.df = data.frame(rbind( IG(ssets.dye[[1]]), IR(ssets.dye[[1]]), II(ssets.dye[[1]]) ))
write.csv(file= file.path(in_dir, 'sesame_mouse_noob_dye.csv'), x=ssets.dye.df, row.names=TRUE)
ssets.betas = lapply(ssets.dye, getBetas)
write.csv(file= file.path(in_dir, 'sesame_mouse_noob_dye_betas.csv'), x=ssets.betas, row.names=TRUE)


#idat_dir = system.file("extdata/", package = "sesameData")
idat_dir = paste0('/Users/mmaxmeister/methylprep/docs/example_data/mouse/')
betas = openSesame(idat_dir, BPPARAM=BiocParallel::MulticoreParam(2))
write.csv(file= file.path(idat_dir, 'open_sesame_mouse_betas.csv'), x=betas, row.names=TRUE)

Rprof(tf <- "log.log", memory.profiling = TRUE)
idat_dir = paste0('/Volumes/LEGX/Barnes/mouse_test/')
betas= openSesame(idat_dir, BPPARAM=BiocParallel::MulticoreParam(2))
Rprof (NULL) ; print(summaryRprof(tf))

system.time( openSesame(idat_dir, BPPARAM=BiocParallel::MulticoreParam(2)) )
