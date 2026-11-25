suppressPackageStartupMessages({library(cmdstanr);library(posterior);library(jsonlite);library(Matrix)})
args <- commandArgs(TRUE); id <- args[1]; attempt <- as.integer(args[2])
root <- normalizePath('outputs/modelos_espectrales'); work <- normalizePath('work/modelos_espectrales')
source(file.path(root,'codigo','funciones_likelihood.R'))
jobs <- fromJSON(file.path(work,'trabajos.json'), simplifyVector=FALSE)
job <- jobs[[which(vapply(jobs,function(j) j$id == id,logical(1)))]]
files <- file.path(job$dir, paste0('intento',attempt,'_cadena',1:4,'.csv'))
read <- read_cmdstan_csv(files,format='draws_array')
draws <- read$post_warmup_draws
diagnostics <- read$post_warmup_sampler_diagnostics
vars <- variables(draws)
core_vars <- vars[!grepl('^v_raw\\[',vars)]
if(job$type != 2) core_vars <- setdiff(core_vars,'expected_K')
core <- subset_draws(draws,variable=core_vars)
summ <- summarise_draws(core)
out <- file.path(root,'resultados',id);dir.create(out,showWarnings=FALSE)
write.csv(as.data.frame(summ),file.path(out,'parametros.csv'),row.names=FALSE)
saveRDS(core,file.path(out,'posterior_parametros.rds'))
dm <- unclass(as_draws_matrix(diagnostics))
energy <- unclass(as_draws_array(subset_draws(diagnostics,variable='energy__')))
bfmi <- vapply(1:4,function(k) mean(diff(energy[,k,1])^2)/var(energy[,k,1]),numeric(1))
max_depth <- if(attempt==1)13 else 14
div <- sum(dm[,'divergent__']);depth <- sum(dm[,'treedepth__']>=max_depth)
rr <- summ$rhat[is.finite(summ$rhat)];ee<-summ$ess_bulk[is.finite(summ$ess_bulk)]
ok <- div == 0 && depth == 0 && max(rr) <= 1.01 && min(ee) >= 400 && min(bfmi) >= .3
diagnostic <- data.frame(id,city=job$city,model_type=job$type,with_covariates=job$with_covariates,
  attempt,chains=4,draws_per_chain=dim(draws)[1],divergences=div,treedepth_hits=depth,
  max_rhat=max(rr),min_bulk_ess=min(ee),min_tail_ess=min(summ$ess_tail,na.rm=TRUE),
  min_ebfmi=min(bfmi),diagnostic_pass=ok)
if(job$type==3) {
  local_summ <- summarise_draws(subset_draws(draws,variable=vars[grepl('^v_raw\\[',vars)]),
                               'rhat','ess_bulk','ess_tail')
  diagnostic$max_local_rhat <- max(local_summ$rhat,na.rm=TRUE)
  diagnostic$min_local_bulk_ess <- min(local_summ$ess_bulk,na.rm=TRUE)
  ok <- ok && diagnostic$max_local_rhat <= 1.01 && diagnostic$min_local_bulk_ess >= 400
  diagnostic$diagnostic_pass <- ok
}
write.csv(diagnostic,file.path(out,'diagnosticos.csv'),row.names=FALSE)
write_json(list(id=id,attempt=attempt,pass=ok,diagnostics=diagnostic),
           file.path(job$dir,'diagnostico_actual.json'),auto_unbox=TRUE,pretty=TRUE)
# Posterior and full fit remain reproducible from the four raw chain files.
saveRDS(list(metadata=read$metadata,diagnostics=diagnostic,chain_files=files),
        file.path(out,'ajuste_metadatos.rds'))
if(!ok && attempt==1) {cat('RETRY NEEDED',id,'\n');quit(status=0)}
d <- readRDS(file.path(root,'datos',paste0(job$city,'_preparado.rds')))
mat <- unclass(as_draws_matrix(draws))
beta_names <- grep('^beta\\[',colnames(mat),value=TRUE)
b_names <- grep('^b\\[',colnames(mat),value=TRUE)
q_names <- grep('^q\\[',colnames(mat),value=TRUE)
v_names <- grep('^v_raw\\[',colnames(mat),value=TRUE)
# Matrix column ordering follows numeric indexes in Stan CSV, not lexical sorting.
X <- if(job$with_covariates)d$X else matrix(numeric(),length(d$y),0)
if(isTRUE(job$free_length))
  X<-cbind(X,log_length_deviation=d$offset-mean(d$offset))
coef_rows <- summ[grepl('^beta\\[',summ$variable),]
if(nrow(coef_rows)) {
  coef_rows$covariate <- colnames(X)[as.integer(sub('beta\\[([0-9]+)\\]','\\1',coef_rows$variable))]
  coef_rows$rate_ratio_median <- exp(coef_rows$median)
  coef_rows$rate_ratio_q05 <- exp(coef_rows$q5);coef_rows$rate_ratio_q95 <- exp(coef_rows$q95)
  write.csv(as.data.frame(coef_rows),file.path(out,'covariables.csv'),row.names=FALSE)
  if(isTRUE(job$free_length)) {
    ls<-coef_rows[coef_rows$covariate=='log_length_deviation',]
    write.csv(data.frame(length_elasticity_median=1+ls$median,
      length_elasticity_q05=1+ls$q5,length_elasticity_q95=1+ls$q95),
      file.path(out,'elasticidad_longitud.csv'),row.names=FALSE)
  }
}
if(length(q_names)) {
  qs <- summarise_draws(subset_draws(draws,variable=q_names),
                        'mean','mcse_mean','ess_bulk','rhat')
  mode_num <- as.integer(sub('q\\[([0-9]+)\\]','\\1',qs$variable))
  pip <- data.frame(mode=mode_num,eigenvalue=d$lambda[mode_num],PIP=qs$mean,
                    MCSE=qs$mcse_mean,ESS=qs$ess_bulk,Rhat=qs$rhat)
  write.csv(pip,file.path(out,'seleccion_modos.csv'),row.names=FALSE)
  kk <- summarise_draws(subset_draws(draws,variable='expected_K'))
  write.csv(as.data.frame(kk),file.path(out,'complejidad.csv'),row.names=FALSE)
  # Exact posterior mixture of Poisson-binomial distributions for the latent K.
  qmat <- mat[,q_names,drop=FALSE]
  pk <- matrix(0,nrow(qmat),ncol(qmat)+1);pk[,1]<-1
  for(j in seq_len(ncol(qmat))) {
    nextpk<-pk*(1-qmat[,j])
    nextpk[,2:ncol(pk)]<-nextpk[,2:ncol(pk),drop=FALSE]+
                         pk[,1:(ncol(pk)-1),drop=FALSE]*qmat[,j]
    pk<-nextpk
  }
  distribution_K <- data.frame(K=0:ncol(qmat),probability=colMeans(pk))
  stopifnot(abs(sum(distribution_K$probability)-1)<1e-8,
    abs(sum(distribution_K$K*distribution_K$probability)-sum(pip$PIP))<1e-8)
  write.csv(distribution_K,file.path(out,'distribucion_K.csv'),row.names=FALSE)
}
# Use every posterior draw for posterior mean structured fields.
if(job$type > 0) {
  weights <- 1/sqrt(1-.9*d$lambda)
  bb <- mat[,b_names,drop=FALSE]
  if(job$type==1) bb <- sweep(bb,2,weights,'*') * mat[,'tau[1]']
  if(job$type==3) bb <- sweep(bb,2,weights,'*') *
    (mat[,'tau[1]']*sqrt(mat[,'rho[1]'])/job$normalizer)
  field_mean <- as.vector(d$B %*% colMeans(bb))
  contribution <- sweep(d$B,2,colMeans(bb),'*')
  write.csv(data.frame(mode=1:20,eigenvalue=d$lambda,mean_coefficient=colMeans(bb),
    rms_contribution=sqrt(colMeans(contribution^2))),file.path(out,'contribuciones_modos.csv'),row.names=FALSE)
} else {field_mean<-rep(0,length(d$y));bb<-matrix(numeric(),nrow(mat),0)}
# Fixed evenly spaced draws in each chain for LOO (up to 2,000 draws total).
per_chain <- min(500L,dim(draws)[1])
ix <- unique(round(seq(1,dim(draws)[1],length.out=per_chain)))
selected <- unclass(as_draws_matrix(subset_draws(draws,iteration=ix)))
S <- nrow(selected);N<-length(d$y)
eta <- matrix(selected[,'alpha'],S,N)
structured_selected <- matrix(0,S,N)
if(ncol(X)) eta <- eta + selected[,beta_names,drop=FALSE] %*% t(X)
if(job$type > 0) {
  bs <- selected[,b_names,drop=FALSE]
  if(job$type==1) bs<-sweep(bs,2,weights,'*')*selected[,'tau[1]']
  if(job$type==3) bs<-sweep(bs,2,weights,'*')*
    (selected[,'tau[1]']*sqrt(selected[,'rho[1]'])/job$normalizer)
  structured_selected <- bs %*% t(d$B)
  eta <- eta + structured_selected
}
eta <- sweep(eta,2,d$offset,'+')
phi <- selected[,'phi']
loglik <- matrix(NA_real_,S,N)
quadrature_nodes_used <- if(job$type==3)30L else 0L
for(i in seq_len(N)) {
  if(job$type != 3) loglik[,i] <- dnbinom(d$y[i],mu=exp(eta[,i]),size=phi,log=TRUE)
  else {
    # Integrate the held-out segment's iid effect against its prior.
    # Do not condition LOO on the fitted observation-specific random effect.
    quad <- if(i==1)statmod::gauss.quad.prob(30,'normal') else quad
    loglik[,i] <- log_nb_quadrature(eta[,i],phi,d$y[i],
       selected[,'tau[1]']*sqrt(1-selected[,'rho[1]']),quad)
  }
}
stopifnot(all(is.finite(loglik)))
if(job$type==3) {
  # Check integration accuracy on extremes and a fixed sample of other segments.
  check_segments <- unique(c(head(order(d$y,decreasing=TRUE),20),
                            round(seq(1,N,length.out=20))))
  quad60 <- statmod::gauss.quad.prob(60,'normal')
  differences <- vapply(check_segments,function(i) {
    ll60<-log_nb_quadrature(eta[,i],phi,d$y[i],
       selected[,'tau[1]']*sqrt(1-selected[,'rho[1]']),quad60)
    max(abs(ll60-loglik[,i]))
  },numeric(1))
  write.csv(data.frame(segment=check_segments,max_abs_loglik_difference_30_vs_60=differences),
            file.path(out,'diagnostico_cuadratura.csv'),row.names=FALSE)
  if(max(differences) > .001) {
    quadrature_nodes_used <- 60L
    cat('Quadrature accuracy: recomputing all integrated log likelihoods with 60 nodes\n')
    for(i in seq_len(N)) {
      loglik[,i]<-log_nb_quadrature(eta[,i],phi,d$y[i],
        selected[,'tau[1]']*sqrt(1-selected[,'rho[1]']),quad60)
    }
  }
}
chain_id <- rep(1:4,each=length(ix))
r_eff <- loo::relative_eff(exp(sweep(loglik,2,apply(loglik,2,max),'-')),chain_id=chain_id,cores=1)
loo_fit <- loo::loo(loglik,r_eff=r_eff,cores=1,save_psis=FALSE)
saveRDS(loo_fit,file.path(out,'loo.rds'))
loo_summary <- data.frame(id,city=job$city,model_type=job$type,with_covariates=job$with_covariates,
  elpd_loo=loo_fit$estimates['elpd_loo','Estimate'],se_elpd=loo_fit$estimates['elpd_loo','SE'],
  p_loo=loo_fit$estimates['p_loo','Estimate'],looic=loo_fit$estimates['looic','Estimate'],
  loo_draws=S,k_threshold=min(1-1/log10(S),.7),
  quadrature_nodes=quadrature_nodes_used,
  pareto_bad=sum(loo_fit$diagnostics$pareto_k>min(1-1/log10(S),.7)),
  pareto_above_1=sum(loo_fit$diagnostics$pareto_k>1),max_pareto_k=max(loo_fit$diagnostics$pareto_k))
write.csv(loo_summary,file.path(out,'prediccion_loo.csv'),row.names=FALSE)
# Posterior predictive checks condition on fitted latent effects, as appropriate for adequacy.
set.seed(job$seed+40000)
take <- unique(round(seq(1,S,length.out=400)))
mu <- exp(eta[take,,drop=FALSE])
if(job$type==3) mu <- mu*exp(selected[take,v_names,drop=FALSE]*
             (selected[take,'tau[1]']*sqrt(1-selected[take,'rho[1]'])))
reps <- matrix(rnbinom(length(mu),mu=as.vector(mu),size=rep(phi[take],N)),nrow=length(take))
statfun <- function(z)c(zero_fraction=mean(z==0),variance=var(z),maximum=max(z),
                        q99=unname(quantile(z,.99)),mean=mean(z))
obs <- statfun(d$y);repstats<-t(apply(reps,1,statfun))
ppc <- data.frame(statistic=names(obs),observed=as.numeric(obs),
  replicated_q05=apply(repstats,2,quantile,.05),replicated_median=apply(repstats,2,median),
  replicated_q95=apply(repstats,2,quantile,.95),bayesian_p_upper=colMeans(sweep(repstats,2,obs,'>=')))
write.csv(ppc,file.path(out,'comprobaciones_predictivas.csv'),row.names=FALSE)
moran <- function(z) {
 z<-z-mean(z);as.numeric(N/sum(d$adjacency)*crossprod(z,d$adjacency%*%z)/sum(z^2))
}
pick<-unique(round(seq(1,nrow(mu),length.out=100)))
Iobs<-Irep<-numeric(length(pick))
for(j in seq_along(pick)) {
 ii<-pick[j];sd<-sqrt(mu[ii,]+mu[ii,]^2/phi[take[ii]])
 Iobs[j]<-moran((d$y-mu[ii,])/sd);Irep[j]<-moran((reps[ii,]-mu[ii,])/sd)
}
write.csv(data.frame(observed_I_median=median(Iobs),replicated_I_q05=quantile(Irep,.05),
 replicated_I_median=median(Irep),replicated_I_q95=quantile(Irep,.95),
 bayesian_p_upper=mean(Irep>=Iobs)),file.path(out,'autocorrelacion_residual.csv'),row.names=FALSE)
write.csv(data.frame(analysis_segment_id=seq_len(N),structured_field_mean=field_mean,
    structured_field_q05=apply(structured_selected,2,quantile,.05),
    structured_field_q95=apply(structured_selected,2,quantile,.95),
    fitted_count_mean=colMeans(mu),fitted_count_q05=apply(mu,2,quantile,.05),
    fitted_count_q95=apply(mu,2,quantile,.95)),file.path(out,'resumen_segmentos.csv'),row.names=FALSE)
saveRDS(list(draw_indices=take,replicated_statistics=repstats,observed_statistics=obs,
             residual_moran_observed=Iobs,residual_moran_replicated=Irep),
        file.path(out,'ppc_resumen.rds'))
cat('COMPLETE',id,'pass=',ok,'elpd=',loo_summary$elpd_loo,'\n')
write_json(list(id=id,complete=TRUE,diagnostic_pass=ok,attempt=attempt),
           file.path(out,'completo.json'),auto_unbox=TRUE,pretty=TRUE)
