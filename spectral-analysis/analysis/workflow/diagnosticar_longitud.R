# Ancillary likelihood diagnostic only; does not replace any Bayesian benchmark.
suppressPackageStartupMessages(library(MASS))
root<-normalizePath('outputs/modelos_espectrales')
res<-list()
for(city in c('Barcelona','Bogota')) {
 d<-readRDS(file.path(root,'datos',paste0(city,'_preparado.rds')))
 dat<-data.frame(y=d$y,d$X,log_length=d$offset)
 fixed<-glm.nb(reformulate(c(colnames(d$X),'offset(log_length)'),response='y'),
               data=dat,control=glm.control(maxit=100))
 free<-glm.nb(reformulate(c(colnames(d$X),'log_length'),response='y'),
              data=dat,control=glm.control(maxit=100))
 est<-coef(summary(free))['log_length',]
 lrt<-2*(as.numeric(logLik(free))-as.numeric(logLik(fixed)))
 res[[city]]<-data.frame(city,length_elasticity=est[1],SE=est[2],
   CI95_low=est[1]-1.96*est[2],CI95_high=est[1]+1.96*est[2],
   coefficient_assumed_by_offset=1,AIC_fixed_offset=AIC(fixed),AIC_free_length=AIC(free),
   likelihood_ratio=lrt,p_value_chisq_1=pchisq(lrt,df=1,lower.tail=FALSE),
   mean_observed=mean(d$y),mean_fitted_fixed=mean(fitted(fixed)),
   mean_fitted_free=mean(fitted(free)),
   note='Frequentist NB with covariates, without spatial field; exploratory offset diagnostic')
 saveRDS(list(fixed_offset=fixed,free_length=free),
         file.path(root,'resultados',paste0(city,'_diagnostico_longitud_glm.rds')))
}
write.csv(do.call(rbind,res),file.path(root,'diagnostico_longitud.csv'),row.names=FALSE)
print(do.call(rbind,res))
