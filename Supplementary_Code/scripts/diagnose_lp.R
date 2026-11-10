suppressPackageStartupMessages({library(cmdstanr);library(posterior)})
p<-commandArgs(TRUE)[1]; f<-as_cmdstan_fit(file.path(p,paste0('chain_',1:4,'.csv')));d<-f$draws('lp__');a<-as.data.frame(summarise_draws(d));write.csv(a,file.path(p,'lp_diagnostics.csv'),row.names=FALSE)
stopifnot(identical(as.integer(dim(as_draws_array(d))[1:2]),c(2000L,4L)))
