suppressPackageStartupMessages({library(cmdstanr);library(posterior)})
target <- commandArgs(TRUE)[1]
files <- file.path(target,paste0('chain_',1:4,'.csv'))
stopifnot(all(file.exists(files)))
fit <- as_cmdstan_fit(files)
d <- fit$diagnostic_summary()
write.csv(data.frame(chain=1:4,divergences=d$num_divergent,treedepth_hits=d$num_max_treedepth,ebfmi=d$ebfmi),file.path(target,'sampler_diagnostics.csv'),row.names=FALSE)
summarize <- function(x) {
 z <- as.data.frame(summarise_draws(x,mean=mean,sd=sd,q025=~unname(quantile(.x,.025)),q975=~unname(quantile(.x,.975)),rhat=rhat,ess_bulk=ess_bulk,ess_tail=ess_tail,mcse_mean=mcse_mean))
 z$mcse_over_sd <- z$mcse_mean/z$sd
 z
}
fixed <- fit$draws(c('alpha','beta','phi','sigma_u'))
write.csv(summarize(fixed),file.path(target,'fixed_diagnostics.csv'),row.names=FALSE)
sigma <- as_draws_array(fit$draws('sigma_u'))[,,1,drop=FALSE]
sigma <- matrix(sigma,nrow=dim(sigma)[1],ncol=dim(sigma)[2])
# Diagnose the actual scaled spatial effects used in the maps, not just raw u.
if (any(grepl('^u\\[',fit$metadata()$model_params))) {
 u <- as_draws_array(fit$draws('u'))
 u <- as_draws_array(sweep(u,c(1,2),sigma,'*'))
} else {
 u <- as_draws_array(fit$draws('spatial_effect'))
}
variables(u) <- sub('u\\[','spatial_effect[',variables(u))
write.csv(summarize(u),file.path(target,'spatial_diagnostics.csv'),row.names=FALSE)
rm(u);gc()
eta <- as_draws_array(fit$draws('eta'))
mu <- exp(eta); rm(eta);gc()
variables(mu) <- sub('eta\\[','mu[',variables(mu))
write.csv(summarize(mu),file.path(target,'fitted_diagnostics.csv'),row.names=FALSE)
rm(mu);gc()
dir.create(file.path(target,'traceplots'),showWarnings=FALSE)
a <- as_draws_array(fixed)
for(v in seq_len(dim(a)[3])) { png(file.path(target,'traceplots',paste0('parameter_',v,'.png')),width=1200,height=700); matplot(seq_len(dim(a)[1]),matrix(as.numeric(a[,,v]),nrow=dim(a)[1],ncol=dim(a)[2]),type='l',lty=1,col=1:4,xlab='Post-warmup iteration',ylab=variables(a)[v],main=variables(a)[v]); dev.off() }
cat('DIAGNOSTICS_COMPLETED',target,'\n')

rows <- lapply(c('fixed','spatial','fitted'),function(kind) {
 z <- read.csv(file.path(target,paste0(kind,'_diagnostics.csv')))
 data.frame(group=kind,n=nrow(z),max_rhat=max(z$rhat),min_ess_bulk=min(z$ess_bulk),min_ess_tail=min(z$ess_tail),max_mcse_over_sd=max(z$mcse_over_sd),n_rhat_ge_1_01=sum(z$rhat>=1.01),n_ess_bulk_lt_400=sum(z$ess_bulk<400),n_ess_tail_lt_400=sum(z$ess_tail<400),n_mcse_over_sd_gt_05=sum(z$mcse_over_sd>.05))
})
write.csv(do.call(rbind,rows),file.path(target,'diagnostic_overview.csv'),row.names=FALSE)
