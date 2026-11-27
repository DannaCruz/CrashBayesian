# Collect every empirical assessment from freshly computed fit outputs.
# Paired predictive uncertainty is recomputed from pointwise LOO vectors, never
# obtained by subtracting marginal standard errors from reference tables.
suppressPackageStartupMessages(library(jsonlite))
a<-commandArgs(TRUE)
if(length(a)!=2)stop('Supply MODEL_RESULTS_ROOT and a NEW OUTPUT_DIRECTORY')
root<-normalizePath(a[1],mustWork=TRUE);out<-a[2]
if(file.exists(out))stop('Refusing an existing metrics output directory')
dir.create(out,recursive=TRUE)
rows<-list();add<-function(finding,fit,metric,value){
 rows[[length(rows)+1]]<<-data.frame(finding=finding,fit_id=fit,metric=metric,value=as.character(value),stringsAsFactors=FALSE)
}
results<-file.path(root,'resultados')
aud<-read.csv(file.path(results,'auditoria_datos.csv'))
for(i in seq_len(nrow(aud)))for(n in names(aud))if(is.numeric(aud[[n]]))add('Table 1 / data audit',aud$city[i],n,aud[[n]][i])
for(city in c('Barcelona','Bogota')){
 spectrum<-read.csv(file.path(results,paste0(city,'_espectro.csv')))
 add('Figure 1 / numerical eigenpair check',city,'maximum_eigen_equation_residual',max(spectrum$residual))
 for(i in seq_len(nrow(spectrum)))add('Figure 1',city,paste0('eigenvalue_',spectrum$index[i]),spectrum$eigenvalue[i])
}
dirs<-list.dirs(results,recursive=FALSE,full.names=TRUE)
dirs<-dirs[file.exists(file.path(dirs,'completo.json'))]
if(length(dirs)!=20)stop('Expected 20 completed empirical summaries')
for(p in dirs){
 id<-basename(p)
 for(file in c('diagnosticos.csv','prediccion_loo.csv','elasticidad_longitud.csv','autocorrelacion_residual.csv')){
  f<-file.path(p,file);if(!file.exists(f))next
  x<-read.csv(f)
  for(n in names(x))if(is.numeric(x[[n]])||is.logical(x[[n]]))add(file,id,n,x[[n]][1])
 }
 pip_file<-file.path(p,'seleccion_modos.csv')
 if(file.exists(pip_file)){
  x<-read.csv(pip_file);pk<-read.csv(file.path(p,'distribucion_K.csv'))
  stopifnot(abs(sum(pk$probability)-1)<1e-8)
  q<-function(prob)pk$K[which(cumsum(pk$probability)>=prob)[1]]
  add('Table 2 / Figure 3',id,'expected_K',sum(x$PIP));add('Table 2',id,'K_q05',q(.05));add('Table 2',id,'K_q95',q(.95))
  add('Table 2 / Figure 3',id,'modes_PIP_above_half',paste(x$mode[x$PIP>.5],collapse=','))
  for(i in seq_len(nrow(x)))for(n in c('PIP','MCSE','ESS','Rhat'))add('Figure 3 / membership diagnostics',id,paste0('mode_',x$mode[i],'_',n),x[[n]][i])
 }
 for(file in c('covariables.csv','comprobaciones_predictivas.csv','contribuciones_modos.csv')){
  f<-file.path(p,file);if(!file.exists(f))next
  x<-read.csv(f);key<-if(file=='covariables.csv')'covariate' else if(file=='comprobaciones_predictivas.csv')'statistic' else 'mode'
  for(i in seq_len(nrow(x)))for(n in names(x))if(is.numeric(x[[n]]))add(file,id,paste(x[[key]][i],n,sep=':'),x[[n]][i])
 }
 seg<-read.csv(file.path(p,'resumen_segmentos.csv'))
 if(grepl('Seleccion_espectral.*longitud_libre$',id)){
  add('Figure 5 / Figure S1',id,'mean_field_RMS',sqrt(mean(seg$structured_field_mean^2)))
  add('Figure S1',id,'segments_interval_above_zero',sum(seg$structured_field_q05>0))
  add('Figure S1',id,'segments_interval_below_zero',sum(seg$structured_field_q95<0))
  add('Figure S1',id,'segments_interval_contains_zero',sum(seg$structured_field_q05<=0&seg$structured_field_q95>=0))
 }
}
contrast<-function(city,left,right,label){
 l<-readRDS(file.path(results,paste(city,left,sep='__'),'loo.rds'))$pointwise[,'elpd_loo']
 r<-readRDS(file.path(results,paste(city,right,sep='__'),'loo.rds'))$pointwise[,'elpd_loo']
 stopifnot(length(l)==length(r),all(is.finite(l)),all(is.finite(r)))
 delta<-l-r
 add('Narrative paired ELPD',paste(city,label,sep=':'),'gain',sum(delta))
 add('Narrative paired ELPD',paste(city,label,sep=':'),'paired_SE',sqrt(length(delta)*var(delta)))
}
for(city in c('Barcelona','Bogota')){
 for(model in c('NB','CAR_espectral','Seleccion_espectral','BYM2_espectral'))contrast(city,paste0(model,'__con_covariables'),paste0(model,'__sin_covariables'),paste0(model,'_attribute_gain'))
 for(model in c('NB','Seleccion_espectral'))contrast(city,paste0(model,'__con_covariables__longitud_libre'),paste0(model,'__con_covariables'),paste0(model,'_length_gain'))
 contrast(city,'Seleccion_espectral__con_covariables__longitud_libre','NB__con_covariables__longitud_libre','mixture_minus_NB_L')
}
all<-do.call(rbind,rows)
write.csv(all,file.path(out,'reported_findings.csv'),row.names=FALSE)
write_json(list(scope='Computed from new fit summaries and pointwise LOO objects; no table values are typed in this collector.',completed_fits=length(dirs),metric_rows=nrow(all)),file.path(out,'metrics_manifest.json'),auto_unbox=TRUE,pretty=TRUE)
cat('Collected',nrow(all),'metric rows from',length(dirs),'fits\n')
