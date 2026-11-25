suppressPackageStartupMessages({library(cmdstanr);library(jsonlite)})
root <- normalizePath('outputs/modelos_espectrales')
work <- normalizePath('work/modelos_espectrales')
jobs <- list(); k <- 0L
# Same segments, offset, basis and priors for the with/without-covariate pairs.
for (type in 0:3) for(city in c('Barcelona','Bogota')) for(cov in c(FALSE,TRUE)) {
  d <- readRDS(file.path(root,'datos',paste0(city,'_preparado.rds')))
  M <- if(type == 0) 0L else 20L
  X <- if(cov) d$X else matrix(numeric(), nrow=length(d$y), ncol=0)
  B <- if(M > 0) d$B else matrix(numeric(), nrow=length(d$y), ncol=0)
  lambda <- if(M > 0) d$lambda else numeric()
  normalizer <- if(M > 0) exp(mean(log(rowSums(sweep(B^2,2,1-0.9*lambda,'/'))))/2) else 1
  dat <- list(N=length(d$y), P=ncol(X), X=X, y=d$y, log_exposure=d$offset,
    M=M, B=B, lambda=lambda, model_type=type, mixture=as.integer(type==2),
    iid=as.integer(type==3), gamma=0.9, spike_sd=0.05, field_normalizer=normalizer)
  label <- c('NB','CAR_espectral','Seleccion_espectral','BYM2_espectral')[type+1]
  id <- paste(city,label,if(cov)'con_covariables' else 'sin_covariables',sep='__')
  dir <- file.path(work,'cadenas',id);dir.create(dir,recursive=TRUE,showWarnings=FALSE)
  data_file <- file.path(dir,'datos.json');write_stan_json(dat,data_file)
  # One deterministic sensible starting point; it is initialization, not an empirical-Bayes prior.
  init <- list(alpha=log((sum(d$y)+0.5)/sum(exp(d$offset))), beta=rep(0,ncol(X)),
               phi=1,b=rep(0,M),tau=if(type>0)array(0.3,dim=1) else numeric(),
               pi_mix=if(type==2)array(0.2,dim=1) else numeric(),
               rho=if(type==3)array(0.5,dim=1) else numeric(),
               v_raw=if(type==3)rep(0,length(d$y)) else numeric())
  init_file <- file.path(dir,'inicial.json');write_stan_json(init,init_file)
  k<-k+1L
  jobs[[k]] <- list(id=id,city=city,type=type,with_covariates=cov,dir=dir,
    data_file=data_file,init_file=init_file,seed=20261120+k,
    normalizer=normalizer,alpha_prior='normal(0,5)',beta_prior='normal(0,0.7)')
}
write_json(jobs,file.path(work,'trabajos.json'),auto_unbox=TRUE,pretty=TRUE)
cat(length(jobs),'jobs ready\n')
