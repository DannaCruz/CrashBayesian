suppressPackageStartupMessages({library(dplyr);library(ggplot2);library(sf)})
root <- normalizePath('outputs/modelos_espectrales')
dirs <- list.dirs(file.path(root,'resultados'),recursive=FALSE,full.names=TRUE)
dirs <- dirs[file.exists(file.path(dirs,'prediccion_loo.csv')) &
               file.exists(file.path(dirs,'resumen_segmentos.csv'))]
extra_dirs<-dirs[grepl('__longitud_libre$',dirs)]
dirs<-dirs[!grepl('__longitud_libre$',dirs)]
read_all <- function(file) bind_rows(lapply(dirs,function(p){
  f<-file.path(p,file);if(!file.exists(f))return(NULL)
  z<-read.csv(f);z$fit_id<-basename(p);z
}))
pred <- read_all('prediccion_loo.csv')
diag <- read_all('diagnosticos.csv')
pred <- pred %>% left_join(diag %>% select(id,diagnostic_pass,max_rhat,min_bulk_ess,
    divergences,min_ebfmi,treedepth_hits),by='id')
model_labels <- c('NB','CAR espectral','Selección espectral','BYM2 espectral')
pred$model <- model_labels[pred$model_type+1]
pred$adjustment <- ifelse(pred$with_covariates,'Con covariables','Sin covariables')
write.csv(pred,file.path(root,'comparacion_modelos.csv'),row.names=FALSE)
write.csv(diag,file.path(root,'diagnosticos_todos.csv'),row.names=FALSE)
extra_diag<-bind_rows(lapply(extra_dirs,function(p) {
 f<-file.path(p,'diagnosticos.csv');if(file.exists(f))read.csv(f) else NULL
}))
write.csv(extra_diag,file.path(root,'diagnosticos_sensibilidad_longitud.csv'),row.names=FALSE)
all_diag<-bind_rows(diag,extra_diag)
all_diag$warnings<-vapply(seq_len(nrow(all_diag)),function(i) {
 z<-all_diag[i,];reasons<-character()
 if(z$divergences>0)reasons<-c(reasons,'divergencias')
 if(z$treedepth_hits>0)reasons<-c(reasons,'profundidad')
 if(z$max_rhat>1.01)reasons<-c(reasons,'Rhat')
 if(z$min_bulk_ess<400)reasons<-c(reasons,'ESS')
 if(z$min_ebfmi<.3)reasons<-c(reasons,'E-BFMI')
 if('max_local_rhat'%in%names(z) && is.finite(z$max_local_rhat) && z$max_local_rhat>1.01)
  reasons<-c(reasons,'Rhat local')
 if('min_local_bulk_ess'%in%names(z) && is.finite(z$min_local_bulk_ess) && z$min_local_bulk_ess<400)
  reasons<-c(reasons,'ESS local')
 if(length(reasons))paste(reasons,collapse=', ') else 'Ninguno'
},character(1))
write.csv(all_diag,file.path(root,'diagnosticos_20_especificaciones.csv'),row.names=FALSE)
coefs <- read_all('covariables.csv')
# Rename the initial internal category label in any older summary snapshot.
if(nrow(coefs))coefs$covariate<-sub('surfacenon_asphalt','surfaceother_recorded',coefs$covariate)
write.csv(coefs,file.path(root,'efectos_covariables.csv'),row.names=FALSE)
pip <- read_all('seleccion_modos.csv')
if(nrow(pip)) {
 pip$city<-sub('__.*','',pip$fit_id)
 pip$adjustment<-ifelse(grepl('__con_covariables$',pip$fit_id),'Con covariables','Sin covariables')
 write.csv(pip,file.path(root,'seleccion_espectral.csv'),row.names=FALSE)
 selection_summary<-bind_rows(lapply(split(pip,pip$fit_id),function(z)
  data.frame(city=z$city[1],adjustment=z$adjustment[1],expected_K=sum(z$PIP),
   modes_PIP_above_05=if(any(z$PIP>.5))paste(z$mode[z$PIP>.5],collapse=', ') else 'Ninguno',
   largest_MCSE=max(z$MCSE))))
 write.csv(selection_summary,file.path(root,'resumen_seleccion.csv'),row.names=FALSE)
 p<-ggplot(pip,aes(mode,PIP,color=adjustment,group=adjustment))+
  geom_hline(yintercept=.5,linetype=2,color='grey60')+geom_line(linewidth=.7)+geom_point(size=1.8)+
  facet_wrap(~city,ncol=1)+scale_x_continuous(breaks=seq(1,20,2))+
  scale_y_continuous(limits=c(0,1),breaks=seq(0,1,.25))+
  scale_color_manual(values=c('Con covariables'='#bc4b27','Sin covariables'='#21618c'))+
  labs(title='Qué modos reciben apoyo después del ajuste',
       subtitle='Offset de longitud fijo; 20 candidatos en la misma base por ciudad',
       x='Índice del modo candidato',y='Probabilidad posterior',color=NULL)+
  theme_minimal(base_size=12)+theme(legend.position='bottom')
 ggsave(file.path(root,'figuras','seleccion_antes_despues.png'),p,width=10,height=7,dpi=180,bg='white')
} else selection_summary<-data.frame()
complexity <- read_all('complejidad.csv')
write.csv(complexity,file.path(root,'complejidad_espectral.csv'),row.names=FALSE)
distK <- read_all('distribucion_K.csv')
write.csv(distK,file.path(root,'distribucion_complejidad_K.csv'),row.names=FALSE)
if(nrow(distK)) {
 k_summary<-bind_rows(lapply(split(distK,distK$fit_id),function(z) {
  z<-z[order(z$K),];cs<-cumsum(z$probability)
  data.frame(fit_id=z$fit_id[1],expected_K=sum(z$K*z$probability),
   K_q05=z$K[which(cs>=.05)[1]],K_median=z$K[which(cs>=.5)[1]],K_q95=z$K[which(cs>=.95)[1]])
 }))
 write.csv(k_summary,file.path(root,'resumen_complejidad_K.csv'),row.names=FALSE)
} else k_summary<-data.frame()
ppc <- read_all('comprobaciones_predictivas.csv')
write.csv(ppc,file.path(root,'comprobaciones_predictivas_todos.csv'),row.names=FALSE)
moran <- read_all('autocorrelacion_residual.csv')
write.csv(moran,file.path(root,'autocorrelacion_residual_todos.csv'),row.names=FALSE)
# Paired ELPD differences use pointwise differences on exactly the same segments.
pairs<-list();kk<-0L
for(city in unique(pred$city))for(type in unique(pred$model_type)) {
 rows<-pred[pred$city==city & pred$model_type==type,]
 if(nrow(rows)!=2 || !all(c(FALSE,TRUE)%in%rows$with_covariates))next
 a<-rows$id[rows$with_covariates];b<-rows$id[!rows$with_covariates]
 la<-readRDS(file.path(root,'resultados',a,'loo.rds'))
 lb<-readRDS(file.path(root,'resultados',b,'loo.rds'))
 diff<-la$pointwise[,'elpd_loo']-lb$pointwise[,'elpd_loo']
 kk<-kk+1L;pairs[[kk]]<-data.frame(city,model=model_labels[type+1],
  elpd_gain_covariates=sum(diff),paired_SE=sqrt(length(diff)*var(diff)),
  diagnostic_pass=all(rows$diagnostic_pass),pareto_bad=sum(rows$pareto_bad))
}
paired<-bind_rows(pairs)
write.csv(paired,file.path(root,'ganancia_covariables.csv'),row.names=FALSE)
length_results<-bind_rows(lapply(extra_dirs,function(p) {
 id<-basename(p);baseid<-sub('__longitud_libre$','',id)
 if(!file.exists(file.path(p,'elasticidad_longitud.csv')))return(NULL)
 ls<-read.csv(file.path(p,'elasticidad_longitud.csv'))
 pr<-read.csv(file.path(p,'prediccion_loo.csv'));dg<-read.csv(file.path(p,'diagnosticos.csv'))
 a<-readRDS(file.path(p,'loo.rds'))
 b<-readRDS(file.path(root,'resultados',baseid,'loo.rds'))
 delta<-a$pointwise[,'elpd_loo']-b$pointwise[,'elpd_loo']
 out<-data.frame(city=pr$city,model=model_labels[pr$model_type+1],ls,
   elpd_gain_free_length=sum(delta),paired_SE=sqrt(length(delta)*var(delta)),
   elpd_loo=pr$elpd_loo,diagnostic_pass=dg$diagnostic_pass,pareto_bad=pr$pareto_bad)
 if(file.exists(file.path(p,'seleccion_modos.csv'))) {
  q<-read.csv(file.path(p,'seleccion_modos.csv'))
  out$expected_K<-sum(q$PIP)
  out$modes_PIP_above_05<-if(any(q$PIP>.5))paste(q$mode[q$PIP>.5],collapse=', ') else 'Ninguno'
 }
 out
}))
write.csv(length_results,file.path(root,'sensibilidad_longitud_bayesiana.csv'),row.names=FALSE)
length_ppc<-bind_rows(lapply(extra_dirs,function(p) {
 f<-file.path(p,'comprobaciones_predictivas.csv');if(!file.exists(f))return(NULL)
 z<-read.csv(f);z$fit_id<-basename(p);z
}))
write.csv(length_ppc,file.path(root,'comprobaciones_sensibilidad_longitud.csv'),row.names=FALSE)
extra_pip<-bind_rows(lapply(extra_dirs,function(p) {
 f<-file.path(p,'seleccion_modos.csv');if(!file.exists(f))return(NULL)
 z<-read.csv(f);z$city<-sub('__.*','',basename(p));z$adjustment<-'Con covariables, longitud libre';z
}))
if(nrow(extra_pip)) {
 original_pip<-pip
 original_pip$adjustment<-paste0(original_pip$adjustment,', offset fijo')
 together<-bind_rows(original_pip,extra_pip)
 write.csv(together,file.path(root,'seleccion_con_sensibilidad_longitud.csv'),row.names=FALSE)
 p<-ggplot(together,aes(mode,PIP,color=adjustment,group=adjustment))+
  geom_hline(yintercept=.5,linetype=2,color='grey60')+geom_line(linewidth=.65)+geom_point(size=1.5)+
  facet_wrap(~city,ncol=1)+scale_y_continuous(limits=c(0,1),breaks=seq(0,1,.25))+
  scale_x_continuous(breaks=seq(1,20,2))+
  scale_color_manual(values=c('Con covariables, longitud libre'='#24784d',
    'Con covariables, offset fijo'='#bc4b27','Sin covariables, offset fijo'='#21618c'))+
  labs(title='Selección espectral y supuesto de longitud',
   subtitle='Misma red y mismos 20 candidatos por ciudad; pertenencia posterior al slab',
   x='Índice del modo candidato',y='Probabilidad posterior',color=NULL)+
  theme_minimal(base_size=12)+theme(legend.position='bottom')+guides(color=guide_legend(nrow=3))
 ggsave(file.path(root,'figuras','seleccion_sensibilidad_longitud.png'),p,
         width=10,height=8,dpi=180,bg='white')
}
if(nrow(paired)) {
 p<-ggplot(paired,aes(elpd_gain_covariates,reorder(model,elpd_gain_covariates)))+
  geom_vline(xintercept=0,color='grey60')+
  geom_errorbar(aes(xmin=elpd_gain_covariates-paired_SE,xmax=elpd_gain_covariates+paired_SE),
                orientation='y',width=.15)+geom_point(size=2.6,color='#21618c')+
  facet_wrap(~city,scales='free_x')+
  labs(title='Cambio predictivo al incorporar covariables',
    subtitle='PSIS-LOO por segmento; barras de un error estándar pareado. No es validación por bloques.',
    x='ELPD con covariables menos ELPD sin covariables',y=NULL)+theme_minimal(base_size=11)
 ggsave(file.path(root,'figuras','ganancia_covariables.png'),p,width=11,height=5,dpi=180,bg='white')
}
for(city in c('Barcelona','Bogota')) {
 d<-readRDS(file.path(root,'datos',paste0(city,'_preparado.rds')))
 spectrum<-read.csv(file.path(root,'resultados',paste0(city,'_espectro.csv')))
 p<-ggplot(spectrum,aes(index,eigenvalue))+geom_line(color='#21618c')+geom_point(size=1.2)+
  geom_vline(xintercept=21.5,linetype=2,color='#bc4b27')+
  labs(title=paste(city,'— extremo superior del espectro'),
       subtitle='Modo 1: Perron omitido. Modos originales 2–21: candidatos del ajuste inicial.',
       x='Índice en el espectro original',y='Autovalor de la adyacencia normalizada')+theme_minimal(base_size=11)
 ggsave(file.path(root,'figuras',paste0(city,'_espectro.png')),p,width=9,height=4.5,dpi=180,bg='white')
 modes<-bind_rows(lapply(1:6,function(j) {
   x<-d$segments;x$mode<-paste('Modo',j);x$raw_eigenvector<-d$U_raw[,j];x
 }))
 p<-ggplot(modes)+geom_sf(aes(color=raw_eigenvector),linewidth=.25)+
  facet_wrap(~mode,ncol=3)+scale_color_gradient2(low='#21618c',mid='grey90',high='#bc4b27',midpoint=0)+
  coord_sf(datum=NA)+theme_void(base_size=11)+
  labs(title=paste(city,'— primeros seis modos candidatos'),
       subtitle='Autovectores originales, antes de centrar y escalar; el signo es arbitrario.',color='Valor')
 ggsave(file.path(root,'figuras',paste0(city,'_modos.png')),p,width=11,height=7,dpi=180,bg='white')
 maps<-list()
 for(cv in c('sin','con')) {
  id<-paste0(city,'__Seleccion_espectral__',cv,'_covariables')
  file<-file.path(root,'resultados',id,'resumen_segmentos.csv')
  if(!file.exists(file))next
  s<-read.csv(file);stopifnot(identical(s$analysis_segment_id,d$segments$analysis_segment_id))
  rr<-d$segments;rr$field<-s$structured_field_mean;rr$adjustment<-paste(cv,'covariables')
  maps[[cv]]<-rr
 }
 if(length(maps)==2) {
  both<-bind_rows(maps)
  p<-ggplot(both)+geom_sf(aes(color=field),linewidth=.35)+facet_wrap(~adjustment,ncol=2)+
   scale_color_gradient2(low='#21618c',mid='grey90',high='#bc4b27',midpoint=0)+
   coord_sf(datum=NA)+theme_void(base_size=11)+
   labs(title=paste(city,'— campo estructurado seleccionado'),
     subtitle='Media posterior del aporte espectral en escala logarítmica; misma escala de color.',color='Campo')
  ggsave(file.path(root,'figuras',paste0(city,'_campo_antes_despues.png')),p,width=11,height=6,dpi=180,bg='white')
 }
 free_field_file<-file.path(root,'resultados',
  paste0(city,'__Seleccion_espectral__con_covariables__longitud_libre'),'resumen_segmentos.csv')
 if(file.exists(free_field_file) && 'con'%in%names(maps)) {
  s<-read.csv(free_field_file)
  rr<-d$segments;rr$field<-s$structured_field_mean;rr$adjustment<-'Longitud estimada'
  fixed<-maps[['con']];fixed$adjustment<-'Offset de longitud fijo'
  p<-ggplot(bind_rows(fixed,rr))+geom_sf(aes(color=field),linewidth=.35)+
   facet_wrap(~adjustment,ncol=2)+
   scale_color_gradient2(low='#21618c',mid='grey90',high='#bc4b27',midpoint=0)+
   coord_sf(datum=NA)+theme_void(base_size=11)+
   labs(title=paste(city,'— sensibilidad del campo espectral a longitud'),
    subtitle='Ambos modelos incluyen las mismas covariables; media posterior del campo estructurado.',color='Campo log')
  ggsave(file.path(root,'figuras',paste0(city,'_campo_sensibilidad_longitud.png')),
         p,width=11,height=6,dpi=180,bg='white')
  rr$sign_support<-ifelse(s$structured_field_q05>0,'Intervalo por encima de 0',
    ifelse(s$structured_field_q95<0,'Intervalo por debajo de 0','Intervalo contiene 0'))
  p<-ggplot(rr)+geom_sf(aes(color=sign_support),linewidth=.35)+
   scale_color_manual(values=c('Intervalo por encima de 0'='#bc4b27',
    'Intervalo por debajo de 0'='#21618c','Intervalo contiene 0'='grey80'))+
   coord_sf(datum=NA)+theme_void(base_size=11)+theme(legend.position='bottom')+
   guides(color=guide_legend(nrow=3))+
   labs(title=paste(city,'— incertidumbre del campo espectral ajustado'),
    subtitle='Covariables y longitud estimada; intervalos marginales del 90%, sin ajuste simultáneo.',color=NULL)
  ggsave(file.path(root,'figuras',paste0(city,'_incertidumbre_campo.png')),
         p,width=8,height=7,dpi=180,bg='white')
 }
 selected_dir<-file.path(root,'resultados',paste0(city,'__Seleccion_espectral__con_covariables'))
 if(file.exists(file.path(selected_dir,'contribuciones_modos.csv'))) {
  cm<-read.csv(file.path(selected_dir,'contribuciones_modos.csv'))
  top<-head(cm$mode[order(cm$rms_contribution,decreasing=TRUE)],6)
  contributions<-bind_rows(lapply(top,function(j) {
   rr<-d$segments;rr$contribution<-d$B[,j]*cm$mean_coefficient[match(j,cm$mode)]
   rr$mode<-factor(paste('Modo',j),levels=paste('Modo',top));rr
  }))
  p<-ggplot(contributions)+geom_sf(aes(color=contribution),linewidth=.3)+facet_wrap(~mode,ncol=3)+
   scale_color_gradient2(low='#21618c',mid='grey90',high='#bc4b27',midpoint=0)+
   coord_sf(datum=NA)+theme_void(base_size=11)+
   labs(title=paste(city,'— mayores contribuciones espectrales con covariables'),
      subtitle='Seis mayores aportes RMS de B_j × E(b_j | datos); no se impone umbral de selección.',
      color='Aporte log')
  ggsave(file.path(root,'figuras',paste0(city,'_contribuciones_ajustadas.png')),
         p,width=11,height=7,dpi=180,bg='white')
 }
}
esc<-function(x) {x<-gsub('&','&amp;',x,fixed=TRUE);x<-gsub('<','&lt;',x,fixed=TRUE);gsub('>','&gt;',x,fixed=TRUE)}
htmltable<-function(x) {
 if(!nrow(x))return('<p>Pendiente.</p>')
 for(n in names(x))if(is.numeric(x[[n]]))x[[n]]<-round(x[[n]],3)
 paste0('<div class="scroll"><table><thead><tr>',paste0('<th>',esc(names(x)),'</th>',collapse=''),
  '</tr></thead><tbody>',paste(apply(x,1,function(row)paste0('<tr>',paste0('<td>',esc(as.character(row)),'</td>',collapse=''),'</tr>')),collapse=''),'</tbody></table></div>')
}
audit<-read.csv(file.path(root,'resultados','auditoria_datos.csv'))
html<-c('<!doctype html><html lang="es"><meta charset="utf-8"><title>Modelos espectrales: corrida inicial</title>',
 '<style>body{font:16px/1.55 system-ui;max-width:1180px;margin:40px auto;padding:0 24px;color:#20313e}h1,h2{line-height:1.2}h2{margin-top:40px}.note{background:#fff3df;padding:18px;border-left:4px solid #bc4b27}table{border-collapse:collapse;font-size:13px;width:100%}th,td{padding:9px;border-bottom:1px solid #ddd;text-align:right}th:first-child,td:first-child{text-align:left}.scroll{overflow:auto}img{max-width:100%;height:auto}a{color:#21618c}</style>',
 '<h1>Modelos espectrales con y sin covariables</h1>',
 paste0('<p>Corrida local iniciada el 20 de noviembre de 2026. Ajustes principales resumidos: ',nrow(pred),' de 16. Sensibilidades a longitud libre: ',nrow(length_results),' de 4.</p>'),
 '<div class="note"><b>Resultado exploratorio.</b> Los conteos se recuperaron de expm1(roads$Y); no se verificaron contra los registros originales. Los atributos OSM no están armonizados con el periodo de los siniestros. PSIS-LOO retiene segmentos vecinos: la validación por bloques y la sensibilidad al número de modos están pendientes. Las ciudades no se comparan en riesgo.</div>',
 '<h2>Datos y soporte</h2>',htmltable(audit[,c('city','original_segments','retained_segments','retained_intersections','retained_count_total','zero_segments','invalid_lanes')]),
 '<p>Se reconstruyó la adyacencia binaria entre segmentos con extremos compartidos, manteniendo vías paralelas. Se seleccionó el componente con más segmentos. En Bogotá se ignoró el vector Y_edges desalineado y se usó roads_ls$Y. La red compartida deja 6.079 segmentos, frente a 6.076 del manuscrito: esta corrida no reproduce exactamente aquella delimitación.</p>',
 '<h2>Especificación</h2>',
 '<p>Binomial negativa para conteos enteros, con offset de longitud en kilómetros. Covariables: tipo de vía (referencia residencial), carriles (referencia dos; uno, tres o más, desconocido) y etiqueta de superficie (referencia asphalt; otra etiqueta, desconocido). Otra etiqueta incluye paved y no significa necesariamente ausencia de asfalto. Un valor de 22 carriles se marcó como desconocido. Las categorías desconocidas son controles auxiliares.</p>',
 '<p>Los modelos espaciales comparten 20 candidatos provenientes de la adyacencia normalizada del grafo de segmentos, omitiendo Perron. Se guardan los autovectores originales para el análisis espectral. Las columnas usadas en los ajustes se centran y escalan por RMS; después de esa transformación no son necesariamente autovectores ni ortogonales. La omisión de Perron puede excluir estructura vinculada al grado.</p>',
 '<p>Priors: intercepto N(0,5), coeficientes de covariables N(0,0.7), dispersión Exp(0.5), escala espacial half-N(0,0.5), gamma fijo 0.9. Selección: pi ~ Beta(1,4), spike N(0,0.05²), slab con desviación sqrt(0.05² + tau²/(1 - gamma lambda)). Este cambio garantiza orden de escalas y requiere distinguir estos resultados de los del manuscrito. BYM2 usa un campo estructurado normalizado a varianza marginal geométrica unitaria, efecto independiente y rho ~ Beta(0.5,0.5).</p>',
 '<p>Cuatro cadenas, 1.000 warm-up y 1.000 muestras por cadena, adapt_delta 0.99 y profundidad 13. Un reajuste automático con 2.000 + 2.000, delta 0.995 y profundidad 14 si hay divergencias, saturación, Rhat &gt; 1.01, ESS bulk &lt; 400 o E-BFMI &lt; 0.3. Las probabilidades de selección se incluyen en los diagnósticos. En BYM2 se revisan también los efectos locales.</p>',
 '<h2>Comparación predictiva y diagnósticos</h2>',htmltable(pred[,c('city','model','adjustment','elpd_loo','se_elpd','p_loo','pareto_bad','max_pareto_k','diagnostic_pass','max_rhat','min_bulk_ess','min_ebfmi','divergences')]),
 paste0('<p>Diagnósticos aprobados: ',sum(all_diag$diagnostic_pass),' de ',nrow(all_diag),' ajustes resumidos, incluyendo sensibilidades. Los avisos de muestreo y los de Pareto corresponden a controles diferentes.</p>'),
 htmltable(all_diag[,c('id','attempt','max_rhat','min_bulk_ess','min_ebfmi','warnings')]),
 '<p>Mayor ELPD indica mejor puntuación predictiva dentro de la misma ciudad. Hasta 2.000 draws para PSIS-LOO, con eficiencia relativa. En BYM2 se integra el efecto independiente del segmento retenido mediante cuadratura normal de 30 puntos, contrastada con 60 en segmentos extremos y otros segmentos de control; se recalcula todo con 60 si la discrepancia supera 0.001 en log densidad. Los casos con Pareto k elevado requieren refit exacto o una validación distinta antes de usar el ranking como definitivo.</p>',
 '<h2>Ganancia al incorporar covariables</h2>',htmltable(paired),
 '<img src="figuras/ganancia_covariables.png" alt="Ganancia predictiva con covariables">',
 '<h2>Selección y complejidad espectral</h2>',
 htmltable(selection_summary),
 '<p>La PIP es pertenencia posterior al componente slab, no probabilidad de coeficiente exactamente distinto de cero. La suma de PIP es la complejidad esperada. El valor 0.5 es una referencia descriptiva; no se usa para truncar las predicciones.</p>',
 '<img src="figuras/seleccion_antes_despues.png" alt="Probabilidades de selección de modos">',htmltable(k_summary),
 '<p>Los intervalos de K provienen de la mezcla posterior exacta de distribuciones Poisson-binomial, no de cuantiles de la suma de probabilidades condicionales.</p>',
 '<h2>Sensibilidad a la proporcionalidad con longitud</h2>',
 '<p>El diagnóstico auxiliar por máxima verosimilitud sugiere que imponer coeficiente uno a log(longitud) puede ser inadecuado. Se añadieron cuatro ajustes bayesianos (NB y selección, ambas ciudades) con las mismas covariables y un coeficiente adicional de log(longitud) centrado. El coeficiente total de longitud tiene prior N(1,0.7). La tabla muestra diferencias pareadas frente al correspondiente modelo con covariables y offset fijo. Son análisis adicionales; no deben mezclarse con la comparación principal de 16 ajustes.</p>',
 htmltable(length_results),
 '<img src="figuras/seleccion_sensibilidad_longitud.png" alt="Selección y sensibilidad a longitud">',
 htmltable(length_ppc),
 '<h2>Comprobaciones predictivas y dependencia residual</h2>',htmltable(ppc),htmltable(moran),
 '<p>Intervalos de referencia del 90% obtenidos con 400 réplicas posteriores. Moran I usa residuos de Pearson centrados y la adyacencia reconstruida; se compara el estadístico observado y replicado en 100 draws. Valores de bayesian_p_upper cercanos a 0 o 1 indican discrepancias que requieren atención. Pasar los diagnósticos del muestreador no implica adecuación del modelo. Son comprobaciones de adecuación condicional, no evidencia causal ni validación fuera de muestra.</p>')
for(city in c('Barcelona','Bogota')) html<-c(html,paste0('<h2>',city,': espectro y mapas</h2>'),
 paste0('<img src="figuras/',city,'_espectro.png" alt="Espectro ',city,'">'),
 paste0('<img src="figuras/',city,'_modos.png" alt="Modos ',city,'">'),
 paste0('<img src="figuras/',city,'_contribuciones_ajustadas.png" alt="Contribuciones ajustadas ',city,'">'),
 paste0('<img src="figuras/',city,'_campo_antes_despues.png" alt="Campos ',city,'">'))
for(city in c('Barcelona','Bogota')) html<-c(html,
 paste0('<img src="figuras/',city,'_campo_sensibilidad_longitud.png" alt="Campo y longitud ',city,'">'),
 paste0('<img src="figuras/',city,'_incertidumbre_campo.png" alt="Incertidumbre ',city,'">'))
html<-c(html,'<h2>Archivos reproducibles</h2><p>comparacion_modelos.csv, ganancia_covariables.csv, seleccion_espectral.csv y efectos_covariables.csv resumen los resultados. datos/ contiene las bases preparadas y los índices de los segmentos originales. codigo/ contiene preparación, modelo Stan, ejecución y resúmenes. Las cadenas completas permanecen en work/modelos_espectrales/cadenas/ junto a sus logs, JSON de datos e inicialización; sus rutas se registran en ajuste_metadatos.rds.</p></html>')
writeLines(html,file.path(root,'informe_resultados.html'))
cat('Report compiled:',nrow(pred),'fits\n')
