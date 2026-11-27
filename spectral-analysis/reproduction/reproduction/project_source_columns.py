"""Project archived public event tables onto the columns needed for reproduction.

No records are deduplicated or reordered here. Every source row is retained with
its row index. Unused personal attributes are omitted; aggregation is a separate,
audited step. Inputs are explicit paths, never an author's personal directory.
"""
from pathlib import Path
import argparse,csv,gzip,hashlib,json,collections

def project(source,output,city):
    source=Path(source);output=Path(output)
    opener=gzip.open if source.suffix=='.gz' else open
    fields=['source_row','event_id','year','date_source','lon','lat']
    years=collections.Counter(); n=0;ids=collections.Counter()
    output.parent.mkdir(parents=True,exist_ok=True)
    # mtime=0 makes the compressed snapshot independent of the export clock.
    with opener(source,'rt',encoding='utf-8-sig',newline='') as src:
        reader=csv.DictReader(src)
        import io
        with output.open('wb') as target, gzip.GzipFile(filename='',mode='wb',fileobj=target,mtime=0) as zipped, io.TextIOWrapper(zipped,encoding='utf-8',newline='') as dest:
            writer=csv.DictWriter(dest,fieldnames=fields,lineterminator='\n');writer.writeheader()
            for n,row in enumerate(reader,1):
                if city=='Barcelona':
                    year=row['NK_Any'];eid=row['Numero_expedient'];lon=row['Longitud_WGS84'];lat=row['Latitud_WGS84']
                    date=f"{int(year):04}-{int(row['Mes_any']):02}-{int(row['Dia_mes']):02}"
                else:
                    year=row['ANO_OCURRENCIA_ACC'];eid=row['FORMULARIO'];lon=row['lon'];lat=row['lat'];date=row['FECHA_OCURRENCIA_ACC']
                writer.writerow(dict(source_row=n,event_id=eid,year=year,date_source=date,lon=lon,lat=lat))
                years[year]+=1;ids[eid]+=1
    return dict(city=city,source_basename=source.name,source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),projection_basename=output.name,projection_sha256=hashlib.sha256(output.read_bytes()).hexdigest(),source_rows=n,unique_identifiers=len(ids),years=dict(sorted(years.items())),columns=fields,projection_policy='All input rows retained in original order; only unused columns removed.')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--barcelona',type=Path,required=True);p.add_argument('--bogota',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.output.exists():p.error('Choose a fresh projection output directory.')
    a.output.mkdir(parents=True)
    meta=[project(a.barcelona,a.output/'Barcelona_records.csv.gz','Barcelona'),project(a.bogota,a.output/'Bogota_records.csv.gz','Bogota')]
    (a.output/'source_projection_manifest.json').write_text(json.dumps(meta,indent=2)+'\n')
    print('Projected source records without filtering, deduplication or row reordering.')
