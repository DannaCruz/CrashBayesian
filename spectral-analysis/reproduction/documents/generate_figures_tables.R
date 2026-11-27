suppressPackageStartupMessages({library(sf);library(ggplot2);library(dplyr);library(grid)})
args<-commandArgs(trailingOnly=TRUE)
# Paths are command-line inputs; no personal source locations are embedded.
if(length(args)!=2) stop('Usage: Rscript generate_figures_tables.R MODEL_RESULTS OUTPUT_DIRECTORY')
src<-normalizePath(args[1],mustWork=TRUE)
dir.create(args[2],recursive=TRUE,showWarnings=FALSE)
out<-normalizePath(args[2],mustWork=TRUE)
dir.create(file.path(out,'figures'),showWarnings=FALSE)
dir.create(file.path(out,'tables'),showWarnings=FALSE)
theme_set(theme_minimal(base_size=11) + theme(text=element_text(size=11), axis.text=element_text(size=11), axis.title=element_text(size=11), strip.text=element_text(size=11), legend.text=element_text(size=11), legend.title=element_text(size=11)))
colors<-c('Unadjusted; fixed length'='#21618c','Road attributes; fixed length'='#bc4b27',
          'Road attributes; estimated length'='#24784d')
savep<-function(name,p,w=6.2,h=5)ggsave(file.path(out,'figures',name),p,width=w,height=h,device=cairo_pdf)
drawpanels<-function(name,plots,nr,nc,w=6.2,h=7) {
 cairo_pdf(file.path(out,'figures',name),width=w,height=h)
 grid.newpage();pushViewport(viewport(layout=grid.layout(nr,nc)))
 for(i in seq_along(plots))print(plots[[i]],vp=viewport(layout.pos.row=(i-1)%/%nc+1,
    layout.pos.col=(i-1)%%nc+1))
 dev.off()
}
cityname<-function(x)ifelse(x=='Bogota','Bogotá',x)
spectra<-bind_rows(lapply(c('Barcelona','Bogota'),function(city){
 x<-read.csv(file.path(src,'resultados',paste0(city,'_espectro.csv')));x$city<-cityname(city);x}))
p<-ggplot(spectra,aes(index,eigenvalue))+geom_line(color='#21618c',linewidth=.55)+
 geom_point(size=1)+geom_vline(xintercept=21.5,linetype=2,color='#bc4b27')+
 facet_wrap(~city,scales='free_y')+
 labs(x='Original spectral index',y='Normalized-adjacency eigenvalue')+
 theme(panel.grid.minor=element_blank())
savep('spectrum.pdf',p,6.2,3.2)
pip<-read.csv(file.path(src,'seleccion_con_sensibilidad_longitud.csv'))
pip$spec<-recode(pip$adjustment,'Sin covariables, offset fijo'='Unadjusted; fixed length',
 'Con covariables, offset fijo'='Road attributes; fixed length',
 'Con covariables, longitud libre'='Road attributes; estimated length')
pip$city<-cityname(pip$city)
p<-ggplot(pip,aes(mode,PIP,color=spec,shape=spec,linetype=spec))+
 geom_hline(yintercept=.5,linetype=2,color='grey55')+geom_line(linewidth=.5)+geom_point(size=1.5)+
 facet_wrap(~city,ncol=1)+scale_color_manual(values=colors)+
 scale_x_continuous(breaks=c(1,5,10,15,20))+
 scale_y_continuous(limits=c(0,1),breaks=c(0,.25,.5,.75,1))+
 labs(x='Candidate mode index',y='Posterior slab-membership probability',color=NULL,shape=NULL,linetype=NULL)+
 theme(legend.position='bottom',legend.text=element_text(size=11),panel.grid.minor=element_blank())+
 guides(color=guide_legend(nrow=3),shape=guide_legend(nrow=3),linetype=guide_legend(nrow=3))
savep('mode_membership.pdf',p,6.2,6.2)
rawplots<-contribplots<-fieldplots<-uncertplots<-list()
for(city in c('Barcelona','Bogota')) {
 d<-readRDS(file.path(src,'datos',paste0(city,'_preparado.rds')))
 for(j in 1:3) {
  x<-d$segments;x$value<-d$U_raw[,j]*sqrt(nrow(x))
  rawplots[[length(rawplots)+1]]<-ggplot(x)+geom_sf(aes(color=value),linewidth=.17)+
   scale_color_gradient2(low='#21618c',mid='grey90',high='#bc4b27',midpoint=0,breaks=scales::breaks_pretty(n=3)) +
   guides(color=guide_colorbar(title.position='top',barwidth=unit(1.7,'in'),barheight=unit(.12,'in'))) +
   coord_sf(datum=NA)+theme_void(base_size=11)+theme(text=element_text(size=11),plot.title=element_text(size=11),legend.text=element_text(size=11),legend.title=element_text(size=11),legend.position='bottom',legend.key.width=unit(.5,'cm'))+
   labs(title=paste0(cityname(city),': mode ',j),color='RMS-scaled value')
 }
 pth<-file.path(src,'resultados',paste0(city,'__Seleccion_espectral__con_covariables__longitud_libre'))
 cm<-read.csv(file.path(pth,'contribuciones_modos.csv'))
 tops<-head(cm$mode[order(cm$rms_contribution,decreasing=TRUE)],3)
 for(j in tops) {
  x<-d$segments;x$value<-d$B[,j]*cm$mean_coefficient[match(j,cm$mode)]
  contribplots[[length(contribplots)+1]]<-ggplot(x)+geom_sf(aes(color=value),linewidth=.17)+
   scale_color_gradient2(low='#21618c',mid='grey90',high='#bc4b27',midpoint=0,breaks=scales::breaks_pretty(n=3)) +
   guides(color=guide_colorbar(title.position='top',barwidth=unit(1.7,'in'),barheight=unit(.12,'in'))) +
   coord_sf(datum=NA)+theme_void(base_size=11)+theme(text=element_text(size=11),plot.title=element_text(size=11),legend.text=element_text(size=11),legend.title=element_text(size=11),legend.position='bottom',legend.key.width=unit(.5,'cm'))+
   labs(title=paste0(cityname(city),': mode ',j),color='Log-count contribution')
 }
 st<-read.csv(file.path(pth,'resumen_segmentos.csv'))
 x<-d$segments;x$value<-st$structured_field_mean
 fieldplots[[city]]<-ggplot(x)+geom_sf(aes(color=value),linewidth=.25)+
  scale_color_gradient2(low='#21618c',mid='grey90',high='#bc4b27',midpoint=0,breaks=scales::breaks_pretty(n=3)) +
   guides(color=guide_colorbar(title.position='top',barwidth=unit(1.7,'in'),barheight=unit(.12,'in'))) +
  coord_sf(datum=NA)+theme_void(base_size=11)+theme(text=element_text(size=11),plot.title=element_text(size=11),legend.text=element_text(size=11),legend.title=element_text(size=11),legend.position='bottom')+
  labs(title=cityname(city),color='Posterior mean log-count contribution')
 x$sign<-ifelse(st$structured_field_q05>0,'Interval above zero',
  ifelse(st$structured_field_q95<0,'Interval below zero','Interval contains zero'))
 uncertplots[[city]]<-ggplot(x)+geom_sf(aes(color=sign),linewidth=.25)+
  scale_color_manual(values=c('Interval above zero'='#bc4b27','Interval below zero'='#21618c','Interval contains zero'='grey80'))+
  coord_sf(datum=NA)+theme_void(base_size=11)+theme(text=element_text(size=11),plot.title=element_text(size=11),legend.text=element_text(size=11),legend.title=element_text(size=11),legend.position='bottom')+
  guides(color=guide_legend(nrow=3))+labs(title=cityname(city),color=NULL)
}
rawplots<-rawplots[c(1,4,2,5,3,6)]
contribplots<-contribplots[c(1,4,2,5,3,6)]
drawpanels('raw_network_modes.pdf',rawplots,3,2,6.2,8.0)
drawpanels('adjusted_contributions.pdf',contribplots,3,2,6.2,8.0)
drawpanels('adjusted_field.pdf',fieldplots,1,2,6.2,3.8)
drawpanels('field_uncertainty.pdf',uncertplots,1,2,6.2,3.8)
# Tables are generated from archived outputs, rather than typed estimates.
fmt<-function(x,d=2)formatC(x,format='f',digits=d)
esc<-function(x)gsub('_','\\_',x,fixed=TRUE)
writebody<-function(name,rows) {
 slash<-intToUtf8(92)
 # Remove the separator before the row terminator, without touching city header cells.
 rows<-vapply(rows,function(row) {
   if (startsWith(row,paste0(slash,'textit'))) return(row)
   while (endsWith(row,strrep(slash,3))) row<-substr(row,1,nchar(row)-1)
   suffix<-paste0(' & ',strrep(slash,2))
   if(endsWith(row,suffix)) row<-paste0(substr(row,1,nchar(row)-nchar(suffix)),' ',strrep(slash,2))
   row<-gsub(paste0(strrep(slash,2),'%'),paste0(slash,'%'),row,fixed=TRUE)
   row
 },character(1))
 writeLines(c(rows,paste0(slash,'bottomrule')),file.path(out,'tables',name))
}
au<-read.csv(file.path(src,'resultados','auditoria_datos.csv'))
writebody('data.tex',c(
 paste('Retained segments',paste(format(au$retained_segments,big.mark=','),collapse=' & '),'\\\\',sep=' & '),
 paste('Retained intersections',paste(format(au$retained_intersections,big.mark=','),collapse=' & '),'\\\\',sep=' & '),
 paste('Recovered total count',paste(format(au$retained_count_total,big.mark=','),collapse=' & '),'\\\\',sep=' & '),
 paste('Represented length (km)',paste(fmt(au$length_km,1),collapse=' & '),'\\\\',sep=' & '),
 paste('Zero-count segments (\\%)',paste(fmt(100*au$zero_segments/au$retained_segments,1),collapse=' & '),'\\\\',sep=' & ')))
rows<-character()
for(city in c('Barcelona','Bogota'))for(spec in c('sin_covariables','con_covariables','con_covariables__longitud_libre')) {
 pth<-file.path(src,'resultados',paste0(city,'__Seleccion_espectral__',spec))
 q<-read.csv(file.path(pth,'seleccion_modos.csv'));pk<-read.csv(file.path(pth,'distribucion_K.csv'))
 cs<-cumsum(pk$probability);lo<-pk$K[which(cs>=.05)[1]];hi<-pk$K[which(cs>=.95)[1]]
 labs<-c(sin_covariables='U',con_covariables='C',con_covariables__longitud_libre='L')
 active<-if(any(q$PIP>.5))paste(q$mode[q$PIP>.5],collapse=', ') else 'None'
 rows<-c(rows,paste(cityname(city),labs[[spec]],fmt(sum(q$PIP)),paste0('[',lo,', ',hi,']'),active,'\\\\',sep=' & '))
}
writebody('membership.tex',rows)
pred<-read.csv(file.path(src,'comparacion_modelos.csv'))
extra<-read.csv(file.path(src,'sensibilidad_longitud_bayesiana.csv'))
rows<-character()
for(city in c('Barcelona','Bogota')) {
 rows<-c(rows,paste0(intToUtf8(92),'textit{',cityname(city),'} & & & & & ',strrep(intToUtf8(92),2)))
 for(model in c('NB','CAR_espectral','Seleccion_espectral','BYM2_espectral'))for(spec in c('sin_covariables','con_covariables','con_covariables__longitud_libre')) {
  id<-paste(city,model,spec,sep='__');pf<-file.path(src,'resultados',id,'prediccion_loo.csv')
  if(!file.exists(pf))next
  pr<-read.csv(pf);dg<-read.csv(file.path(dirname(pf),'diagnosticos.csv'))
  lab<-c(NB='NB',CAR_espectral='Dense spectral',Seleccion_espectral='Mixture spectral',BYM2_espectral='Spectral BYM2')[[model]]
  specshort<-c(sin_covariables='U',con_covariables='C',con_covariables__longitud_libre='L')[[spec]]
  rows<-c(rows,paste(lab,specshort,fmt(pr$elpd_loo,1),fmt(pr$se_elpd,1),pr$pareto_bad,
      if(dg$diagnostic_pass)'Pass' else 'Flagged','\\\\',sep=' & '))
 }
 rows<-c(rows,'\\addlinespace')
}
writebody('prediction.tex',rows)
names_cov<-c('highwaytertiary','highwaysecondary','highwayprimary','lanesone','lanesthree_plus','surfaceother_recorded')
labels_cov<-c('Tertiary vs residential road','Secondary vs residential road','Primary vs residential road',
 'One vs two lanes','Three or more vs two lanes','Other surface tag vs asphalt')
effects<-lapply(c('Barcelona','Bogota'),function(city)read.csv(file.path(src,'resultados',
 paste0(city,'__Seleccion_espectral__con_covariables__longitud_libre'),'covariables.csv')))
rows<-vapply(seq_along(names_cov),function(i) {
 vals<-vapply(effects,function(z){r<-z[z$covariate==names_cov[i],];paste0(fmt(r$rate_ratio_median),' [',fmt(r$rate_ratio_q05),', ',fmt(r$rate_ratio_q95),']')},character(1))
 paste(labels_cov[i],vals[1],vals[2],'\\\\',sep=' & ')
},character(1))
writebody('effects.tex',rows)
cat('English figures and source-linked tables generated\n')
