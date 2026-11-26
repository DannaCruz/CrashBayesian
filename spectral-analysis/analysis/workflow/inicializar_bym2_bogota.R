# Starting points only; likelihood and priors remain those in the delivered Stan source.
suppressPackageStartupMessages({library(cmdstanr);library(jsonlite)})
root<-normalizePath('outputs/modelos_espectrales');work<-normalizePath('work/modelos_espectrales')
jobs<-fromJSON(file.path(work,'trabajos.json'),simplifyVector=FALSE)
for(j in jobs)if(j$type==3 && j$city=='Bogota') {
 d<-readRDS(file.path(root,'datos','Bogota_preparado.rds'))
 refid<-sub('BYM2_espectral','CAR_espectral',j$id,fixed=TRUE)
 ref<-read.csv(file.path(root,'resultados',refid,'parametros.csv'))
 val<-function(n)ref$median[match(n,ref$variable)]
 tau<-1.2;rho<-.1;sigma_u<-tau*sqrt(1-rho)
 alpha<-val('alpha')-.5*sigma_u^2
 beta<-if(j$with_covariates)val(paste0('beta[',1:ncol(d$X),']')) else numeric()
 b<-val(paste0('b[',1:20,']'))*val('tau[1]')*j$normalizer/(tau*sqrt(rho))
 eta<-d$offset+alpha+as.vector(d$B%*%(b/sqrt(1-.9*d$lambda)))*tau*sqrt(rho)/j$normalizer
 if(length(beta))eta<-eta+as.vector(d$X%*%beta)
 v<-pmax(-3,pmin(3,.5*(log(d$y+.5)-eta)/sigma_u))
 old<-file.path(j$dir,'inicial_original.json')
 if(!file.exists(old))file.copy(j$init_file,old)
 write_stan_json(list(alpha=alpha,beta=beta,phi=3,b=b,tau=array(tau,1),
    rho=array(rho,1),pi_mix=numeric(),v_raw=v),j$init_file)
 cat('Updated initialization for',j$id,'\n')
}
