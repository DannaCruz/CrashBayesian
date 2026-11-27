# Reconstruct segment responses from frozen source records and road snapshots.
# Usage: Rscript rebuild_from_sources.R COUNT_UNIT SOURCE_DATA OUTPUT_DIRECTORY
# COUNT_UNIT is historical_records (faithful original reconstruction) or
# unique_crashes (one event per year/incident identifier on the analysis support).
suppressPackageStartupMessages({library(sf);library(sfnetworks);library(igraph);library(jsonlite)})
args<-commandArgs(TRUE)
if(length(args)!=3 || !args[1]%in%c('historical_records','unique_crashes'))
 stop('Supply historical_records or unique_crashes, SOURCE_DATA and a NEW OUTPUT_DIRECTORY')
mode<-args[1];src<-normalizePath(args[2],mustWork=TRUE);dest<-args[3]
if(dir.exists(dest)||file.exists(dest))stop('Refusing existing output directory')
dir.create(dest,recursive=TRUE);dest<-normalizePath(dest)
dir.create(file.path(dest,'audits'));dir.create(file.path(dest,'networks'))
audits<-list()
for(city in c('Barcelona','Bogota')) {
 raw<-read.csv(gzfile(file.path(src,paste0(city,'_records.csv.gz'))),
  colClasses=c('integer','character','integer','character','numeric','numeric'),check.names=FALSE)
 n_source<-nrow(raw)
 stopifnot(identical(names(raw),c('source_row','event_id','year','date_source','lon','lat')),
  all(!is.na(raw$event_id)&nzchar(raw$event_id)),all(is.finite(raw$year)))
 valid<-is.finite(raw$lon)&is.finite(raw$lat)&abs(raw$lon)<=180&abs(raw$lat)<=90
 invalid_rows<-sum(!valid);raw<-raw[valid,]
 key<-paste(raw$year,raw$event_id,sep=':')
 # Do not choose a location arbitrarily when duplicated incident IDs disagree.
 coord<-paste(format(raw$lon,digits=15),format(raw$lat,digits=15),sep=':')
 first<-match(key,key);conflicts<-sum(coord!=coord[first])
 if(conflicts)stop(paste(city,'has conflicting locations for an incident identifier; inspect source before analysis'))
 duplicate_records<-sum(duplicated(key))
 if(mode=='unique_crashes')raw<-raw[!duplicated(key),]
 roads<-readRDS(file.path(src,'networks',paste0(city,'_analysis_roads.rds')))
 stopifnot(inherits(roads,'sf'),!anyDuplicated(roads$osm_id),!('Y'%in%names(roads)))
 pts<-st_as_sf(raw,coords=c('lon','lat'),crs=4326,remove=FALSE)
 pts<-st_transform(pts,st_crs(roads))
 if(mode=='historical_records' && city=='Bogota') {
  # The archived workflow assigned all source rows to the complete source road
  # snapshot, aggregated by OSM ID, then retained selected road IDs. Preserve this
  # path only as a labelled reconstruction of the original record-count analysis.
  assignment_roads<-readRDS(file.path(src,'networks','Bogota_assignment_roads.rds'))
  domain<-st_as_sfc(st_bbox(assignment_roads))
 } else {
  assignment_roads<-roads
  # Unique-crash inference uses the actual analysis-support domain before nearest
  # assignment. It must not import counts from clipped-away portions of OSM roads.
  if(city=='Bogota')domain<-st_transform(readRDS(file.path(src,'networks','Bogota_domain.rds')),st_crs(roads))
  else domain<-st_as_sfc(st_bbox(roads))
 }
 inside<-lengths(st_intersects(pts,domain))>0
 outside_rows<-sum(!inside);pts<-pts[inside,];raw<-raw[inside,]
 nearest<-st_nearest_feature(pts,assignment_roads)
 assigned_osm<-as.character(assignment_roads$osm_id[nearest])
 target<-match(assigned_osm,as.character(roads$osm_id))
 counts<-tabulate(target,nbins=nrow(roads))
 stopifnot(all(counts>=0),sum(counts)==sum(!is.na(target)))
 distance<-as.numeric(st_distance(pts,assignment_roads[nearest,],by_element=TRUE))
 # Fixed geometry and input order define graph support; no response is consulted.
 roads$Y<-log1p(counts)
 net<-as_sfnetwork(roads,directed=FALSE);g<-as.igraph(net)
 stopifnot(ecount(g)==nrow(roads),identical(st_as_binary(st_geometry(roads)),
  st_as_binary(st_geometry(st_as_sf(net,'edges')))))
 Y_edges<-counts
 response_provenance<-list(count_unit=mode,city=city,source_rows=n_source,
  invalid_coordinates=invalid_rows,duplicate_valid_rows=duplicate_records,
  year_min=min(raw$year),year_max=max(raw$year),
  response_definition=if(mode=='unique_crashes')'One source incident identifier per year, assigned once within the stated analysis domain' else if(city=='Bogota')'Person-involvement source rows, preserving the archived aggregation' else 'Accident-report source rows, preserving the archived aggregation')
 if(city=='Bogota') {
  roads_ls<-roads;save(roads_ls,net,g,Y_edges,response_provenance,file=file.path(dest,'networks','Bogota_network_ready_for_models.RData'))
 } else {
  roads_sf<-roads;save(roads_sf,net,g,Y_edges,response_provenance,file=file.path(dest,'networks','Eixample_network_ready_for_models.RData'))
 }
 assignments<-data.frame(source_row=raw$source_row,event_id=raw$event_id,year=raw$year,
  source_road_row=nearest,source_osm_id=assigned_osm,analysis_source_row=target,
  assignment_distance_m=distance)
 con<-gzfile(file.path(dest,'audits',paste0(city,'_assignments.csv.gz')),'wt')
 write.csv(assignments,con,row.names=FALSE);close(con)
 audit<-data.frame(city,count_unit=mode,source_rows=n_source,invalid_coordinates=invalid_rows,
  duplicate_valid_rows=duplicate_records,rows_excluded_by_domain=outside_rows,
  input_rows_assigned=nrow(raw),rows_on_selected_road_ids=sum(!is.na(target)),
  source_road_features=nrow(assignment_roads),analysis_source_segments=nrow(roads),
  source_support_count_total=sum(counts),source_year_min=min(raw$year),source_year_max=max(raw$year),
  median_assignment_distance_m=median(distance),max_assignment_distance_m=max(distance))
 audits[[city]]<-audit
 write_json(response_provenance,file.path(dest,'audits',paste0(city,'_provenance.json')),auto_unbox=TRUE,pretty=TRUE)
 cat(city,mode,':',sum(counts),'counts on',nrow(roads),'source segments\n')
}
write.csv(do.call(rbind,audits),file.path(dest,'audits','source_workflow_audit.csv'),row.names=FALSE)
