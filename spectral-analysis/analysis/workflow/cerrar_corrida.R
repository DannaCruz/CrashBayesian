suppressPackageStartupMessages({library(cmdstanr);library(posterior);library(jsonlite)})
root<-normalizePath('outputs/modelos_espectrales')
jobs<-fromJSON('work/modelos_espectrales/trabajos.json',simplifyVector=FALSE)
checks<-list()
for(j in jobs) {
 p<-file.path(root,'resultados',j$id)
 stopifnot(file.exists(file.path(p,'completo.json')))
 meta<-readRDS(file.path(p,'ajuste_metadatos.rds'))
 if(!file.exists(file.path(p,'posterior_parametros.rds'))) {
  z<-read_cmdstan_csv(meta$chain_files,format='draws_array')
  core_vars<-read.csv(file.path(p,'parametros.csv'))$variable
  saveRDS(subset_draws(z$post_warmup_draws,variable=core_vars),
          file.path(p,'posterior_parametros.rds'))
 }
 d<-readRDS(file.path(root,'datos',paste0(j$city,'_preparado.rds')))
 seg<-read.csv(file.path(p,'resumen_segmentos.csv'))
 stopifnot(nrow(seg)==length(d$y),
  identical(seg$analysis_segment_id,d$segments$analysis_segment_id),
  all(is.finite(seg$fitted_count_mean)),all(seg$fitted_count_mean>0))
 loo<-readRDS(file.path(p,'loo.rds'));ls<-read.csv(file.path(p,'prediccion_loo.csv'))
 stopifnot(nrow(loo$pointwise)==length(d$y),
  abs(sum(loo$pointwise[,'elpd_loo'])-ls$elpd_loo)<1e-6)
 if(j$type==2) {
  pip<-read.csv(file.path(p,'seleccion_modos.csv'))
  pk<-read.csv(file.path(p,'distribucion_K.csv'))
  stopifnot(nrow(pip)==20,all(pip$PIP>=0 & pip$PIP<=1),
    all(is.finite(pip$MCSE)),abs(sum(pk$probability)-1)<1e-8,
    abs(sum(pk$K*pk$probability)-sum(pip$PIP))<1e-8)
 }
 checks[[j$id]]<-data.frame(id=j$id,verified_segments=nrow(seg),
    chains=length(meta$chain_files),verified_loo=TRUE,verified_selection=if(j$type==2)TRUE else NA,
    sampler_diagnostic_pass=meta$diagnostics$diagnostic_pass)
}
write.csv(do.call(rbind,checks),file.path(root,'verificacion_resultados.csv'),row.names=FALSE)
pkgs<-c('sf','sfnetworks','igraph','Matrix','RSpectra','cmdstanr','posterior','loo','statmod','matrixStats','jsonlite','dplyr','ggplot2')
versions<-data.frame(package=pkgs,version=vapply(pkgs,function(p)as.character(packageVersion(p)),character(1)))
write.csv(versions,file.path(root,'versiones_paquetes.csv'),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(root,'sesion_R.txt'))
files<-list.files(file.path(root,'codigo'),full.names=TRUE,pattern='\\.(R|py|stan)$')
manifest<-data.frame(file=basename(files),md5=unname(tools::md5sum(files)))
write.csv(manifest,file.path(root,'huellas_codigo.csv'),row.names=FALSE)
cat('Verified',length(jobs),'completed fits\n')
