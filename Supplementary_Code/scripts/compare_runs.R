compare <- function(first,second,prefix) {
 overview <- list()
 for(kind in c('fixed','spatial','fitted')) {
  a <- read.csv(file.path(first,paste0(kind,'_diagnostics.csv')))
  b <- read.csv(file.path(second,paste0(kind,'_diagnostics.csv')))
  stopifnot(!anyDuplicated(a$variable),!anyDuplicated(b$variable),setequal(a$variable,b$variable))
  z <- merge(a,b,by='variable',suffixes=c('_run1','_run2'))
  z$mean_difference <- z$mean_run2-z$mean_run1
  z$combined_mcse <- sqrt(z$mcse_mean_run1^2+z$mcse_mean_run2^2)
  z$difference_over_mcse <- z$mean_difference/z$combined_mcse
  z$difference_over_sd <- z$mean_difference/sqrt((z$sd_run1^2+z$sd_run2^2)/2)
  write.csv(z,paste0(prefix,'_',kind,'_stability.csv'),row.names=FALSE)
  top_n<-max(1,ceiling(.1*nrow(z)))
  overview[[kind]]<-data.frame(group=kind,n=nrow(z),max_abs_difference_over_sd=max(abs(z$difference_over_sd)),median_abs_difference_over_sd=median(abs(z$difference_over_sd)),n_abs_difference_over_mcse_gt_3=sum(abs(z$difference_over_mcse)>3),mean_rank_correlation=if(kind=='fixed') NA else cor(z$mean_run1,z$mean_run2,method='spearman'),top_decile_overlap=if(kind=='fixed') NA else length(intersect(order(z$mean_run1,decreasing=TRUE)[1:top_n],order(z$mean_run2,decreasing=TRUE)[1:top_n]))/top_n)
 }
 write.csv(do.call(rbind,overview),paste0(prefix,'_stability_overview.csv'),row.names=FALSE)
}
args<-commandArgs(TRUE)
if(length(args)==3) {
 compare(args[1],args[2],args[3])
} else if(length(args)==0) {
 for(s in c('intersection','segment')) compare(file.path('runs',paste0(s,'_1')),file.path('runs',paste0(s,'_2')),s)
} else stop('Usage: compare_runs.R first_run_directory second_run_directory output_prefix')
