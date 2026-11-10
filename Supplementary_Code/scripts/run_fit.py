#!/usr/bin/env python3
"""Validate inputs or explicitly reproduce one fit. Never overwrites a run."""
import argparse,json,os,pathlib,subprocess,math
root=pathlib.Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');p.add_argument('--label');p.add_argument('--output');a=p.parse_args()
runs=json.loads((root/'runs.json').read_text())
for r in runs:
 d=json.loads((root/r['data']).read_text());assert len(d['Y'])==d['N']==len(d['X'])==len(d['log_offset'])
 assert all(len(x)==d['K'] for x in d['X']);assert all(y>=0 and int(y)==y for y in d['Y'])
 assert all(math.isfinite(v) for v in d['log_offset']);assert (root/r['model']).is_file()
 if r['init'] not in ['0','0.1']:
  for ch in range(1,5):json.loads((root/r['init']/f'chain_{ch}.json').read_text())
 print(r['label'],d['N'],d['K'],sum(d['Y']),[r['base_seed']+c for c in range(1,5)])
if a.check:raise SystemExit(0)
if not a.label or not a.output:p.error('Use --check OR --label NAME --output NEW_DIRECTORY')
r=next((r for r in runs if r['label']==a.label),None)
if r is None:p.error('Unknown label')
target=pathlib.Path(a.output).resolve()
if target.exists():raise SystemExit('Refusing to overwrite existing output')
cmdstan=os.environ.get('CMDSTAN')
if not cmdstan:raise SystemExit('Set CMDSTAN to your installed CmdStan directory')
exe=(root/r['model']).with_suffix('')
subprocess.run(['make','-C',cmdstan,str(exe)],check=True)
target.mkdir(parents=True)
for ch in range(1,5):
 init=r['init'] if r['init'] in ['0','0.1'] else str(root/r['init']/f'chain_{ch}.json')
 cmd=[str(exe),'sample',f"num_warmup={r['warmup']}",f"num_samples={r['samples']}",'adapt',f"delta={r['adapt_delta']}",'algorithm=hmc','engine=nuts',f"max_depth={r['max_depth']}",f"stepsize={r['step_size']}",'data','file='+str(root/r['data']),'init='+init,'random',f"seed={r['base_seed']+ch}",f'id={ch}','output','file='+str(target/f'chain_{ch}.csv'),f"sig_figs={r['sig_figs']}",'refresh=100']
 with (target/f'chain_{ch}.log').open('w') as f:subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,check=True)
if 'nonspatial' not in r['label']:
 for name in ['diagnose.R','diagnose_lp.R']:subprocess.run(['Rscript',str(root/'scripts'/name),str(target)],check=True)
print('SAMPLING_COMPLETE; assess diagnostics before interpreting the results.')
