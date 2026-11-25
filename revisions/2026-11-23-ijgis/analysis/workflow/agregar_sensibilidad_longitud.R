suppressPackageStartupMessages({library(cmdstanr);library(jsonlite)})
root<-normalizePath('outputs/modelos_espectrales');work<-normalizePath('work/modelos_espectrales')
jobs<-fromJSON(file.path(work,'trabajos.json'),simplifyVector=FALSE)
jobs<-jobs[!vapply(jobs,function(j)isTRUE(j$free_length),logical(1))]
extras<-list();i<-0L
for(type in c(0L,2L))for(city in c('Barcelona','Bogota')) {
 baseid<-paste(city,c('NB','CAR_espectral','Seleccion_espectral','BYM2_espectral')[type+1],
               'con_covariables',sep='__')
 j<-jobs[[which(vapply(jobs,function(z)z$id==baseid,logical(1)))]]
 dat<-fromJSON(j$data_file)
 d<-readRDS(file.path(root,'datos',paste0(city,'_preparado.rds')))
 dat$X<-cbind(d$X,log_length_deviation=d$offset-mean(d$offset));dat$P<-ncol(dat$X)
 j$id<-paste0(j$id,'__longitud_libre');j$free_length<-TRUE
 j$dir<-file.path(work,'cadenas',j$id);dir.create(j$dir,recursive=TRUE,showWarnings=FALSE)
 j$data_file<-file.path(j$dir,'datos.json')
 M<-if(type==0)0L else 20L
 dat$B<-if(M>0)d$B else matrix(numeric(),length(d$y),0)
 dat$lambda<-if(M>0)d$lambda else numeric()
 write_stan_json(dat,j$data_file)
 init<-list(alpha=log((sum(d$y)+.5)/sum(exp(d$offset))),
    beta=c(rep(0,ncol(d$X)),-.4),phi=1,b=rep(0,M),
    tau=if(type>0)array(.3,1) else numeric(),
    pi_mix=if(type==2)array(.2,1) else numeric(),rho=numeric(),v_raw=numeric())
 j$init_file<-file.path(j$dir,'inicial.json');write_stan_json(init,j$init_file)
 i<-i+1L;j$seed<-20261200+i;extras[[i]]<-j
}
write_json(c(jobs,extras),file.path(work,'trabajos.json'),auto_unbox=TRUE,pretty=TRUE)
cat(length(jobs),'main fits plus',length(extras),'length-sensitivity fits\n')
