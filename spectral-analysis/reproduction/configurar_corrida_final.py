"""Validate and apply the archived settings without running a sampler."""
from pathlib import Path
import csv,json,math,shutil

def validate_reference_settings(model_root):
    model_root=Path(model_root)
    runs=json.loads((model_root/'configuraciones_finales.json').read_text())
    assert len(runs)==20 and len({r['id'] for r in runs})==20
    sizes={}
    for city in ['Barcelona','Bogota']:
        with (model_root/'datos'/f'{city}_segmentos.csv').open() as f:
            sizes[city]=sum(1 for _ in csv.DictReader(f))
    def finite(x):
        if isinstance(x,list):return all(finite(v) for v in x)
        if isinstance(x,dict):return all(finite(v) for v in x.values())
        return isinstance(x,(int,float)) and math.isfinite(x)
    for r in runs:
        name=Path(r['initialization']).name
        init=json.loads((model_root/'inicializaciones_finales'/name).read_text())
        assert finite(init),r['id']
        assert len(init['beta'])==(8 if r['with_covariates'] else 0)+int(r['free_length']),r['id']
        assert len(init['b'])==(20 if r['model_type'] else 0),r['id']
        assert len(init['v_raw'])==(sizes[r['city']] if r['model_type']==3 else 0),r['id']
        assert r['chains']==4 and r['sampling'] in [1000,2000]
        assert r['warmup']==r['sampling']
        assert r['max_depth']==(14 if r['sampling']==2000 else 13)
        assert r['adapt_delta']==(.995 if r['sampling']==2000 else .99)
    return runs

def configure(model_root,jobs_file):
    model_root=Path(model_root);jobs_file=Path(jobs_file)
    runs={r['id']:r for r in validate_reference_settings(model_root)}
    jobs=json.loads(jobs_file.read_text())
    assert {j['id'] for j in jobs}==set(runs)
    for j in jobs:
        r=runs[j['id']]
        assert j['type']==r['model_type'] and j['seed']==r['seed'],j['id']
        assert j['with_covariates']==r['with_covariates'],j['id']
        assert bool(j.get('free_length',False))==r['free_length'],j['id']
        shutil.copy2(model_root/'inicializaciones_finales'/Path(r['initialization']).name,j['init_file'])
        j['final_settings']={k:r[k] for k in ['warmup','sampling','adapt_delta','max_depth']}
    jobs_file.write_text(json.dumps(jobs,indent=2))
    return len(jobs)
