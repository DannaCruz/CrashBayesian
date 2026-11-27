"""Compare computed numerical summaries with the deposited reference fits."""
import argparse,csv,json,math
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--fits',type=Path,required=True);p.add_argument('--reference',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
checks=[]
for ref in sorted(a.reference.glob('*/*.csv')):
 if not (ref.parent/'completo.json').exists():continue
 fresh=a.fits/ref.relative_to(a.reference)
 if not fresh.is_file():raise SystemExit('Missing computed summary '+str(fresh))
 with ref.open() as f:r=list(csv.DictReader(f))
 with fresh.open() as f:n=list(csv.DictReader(f))
 if len(r)!=len(n):raise SystemExit('Row count mismatch '+str(fresh))
 numeric=0;maxdiff=0.;failed=0
 for x,y in zip(r,n):
  if not x.keys()<=y.keys():raise SystemExit('Missing reference columns '+str(fresh))
  for key in x:
   try:v=float(x[key]);w=float(y[key])
   except (ValueError,TypeError):continue
   if math.isnan(v) and math.isnan(w):continue
   numeric+=1
   if v==w:continue
   delta=abs(v-w);maxdiff=max(maxdiff,delta)
   if not math.isclose(v,w,rel_tol=1e-7,abs_tol=1e-8):failed+=1
 checks.append(dict(file=str(ref.relative_to(a.reference)),numeric_cells=numeric,max_absolute_difference=maxdiff,cells_outside_tolerance=failed))
report=dict(computed_fit_count=len(list(a.fits.glob('*/completo.json'))),compared_files=len(checks),numeric_cells=sum(x['numeric_cells'] for x in checks),cells_outside_tolerance=sum(x['cells_outside_tolerance'] for x in checks),checks=checks)
a.output.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k!='checks'},indent=2))
if report['computed_fit_count']!=20:raise SystemExit('Expected 20 fits')
