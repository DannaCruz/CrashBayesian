"""Run the empirical fits via CmdStan CLI; all subprocesses and files stay local.

Run from the workspace root after preparar_datos.R, compilation and crear_trabajos.R.
Each fit uses four chains, 1,000 warm-up and 1,000 sampling iterations per chain.
One diagnostic-triggered retry uses 2,000 + 2,000 iterations and target acceptance .995.
"""
from pathlib import Path
import concurrent.futures
import argparse
import json
import subprocess
import shutil
import time

ROOT = Path.cwd()
WORK = ROOT / 'work/modelos_espectrales'
OUT = ROOT / 'modelos_espectrales'
EXE = WORK / 'stan/modelo_nb_espectral'
parser = argparse.ArgumentParser()
parser.add_argument('--only-length',action='store_true')
parser.add_argument('--final-settings',action='store_true',help='Use archived final initialization and sampling settings; do not retry with another budget.')
parser.add_argument('--workers',type=int,default=4)
args = parser.parse_args()
JOBS = json.loads((WORK / 'trabajos.json').read_text())
if args.only_length:
    JOBS = [j for j in JOBS if j.get('free_length',False)]
STATUS = WORK / ('estado_longitud.json' if args.only_length else 'estado.json')

def status(**kw):
    kw['updated_at'] = time.strftime('%Y-%m-%dT%H:%M:%S%z')
    STATUS.write_text(json.dumps(kw, indent=2))

def chain(job, attempt, chain_id):
    out = Path(job['dir'])
    csv = out / f'intento{attempt}_cadena{chain_id}.csv'
    log = out / f'intento{attempt}_cadena{chain_id}.log'
    final = job.get('final_settings') if args.final_settings else None
    n = final['sampling'] if final else (1000 if attempt == 1 else 2000)
    warmup = final['warmup'] if final else n
    target = str(final['adapt_delta']) if final else ('.99' if attempt == 1 else '.995')
    depth = str(final['max_depth']) if final else ('13' if attempt == 1 else '14')
    cmd = [str(EXE), 'sample', f'num_samples={n}', f'num_warmup={warmup}',
           'save_warmup=0', 'thin=1', 'adapt', f'delta={target}',
           'algorithm=hmc', 'engine=nuts', f'max_depth={depth}',
           'data', f"file={job['data_file']}",
           f"init={job['init_file']}", 'random', f"seed={job['seed']}",
           f'id={chain_id}', 'output', f'file={csv}', 'refresh=200', 'sig_figs=8']
    with log.open('w') as stream:
        result = subprocess.run(cmd, stdout=stream, stderr=subprocess.STDOUT)
    if result.returncode:
        raise RuntimeError(f'Chain {chain_id} failed for {job["id"]}: {log.read_text()[-3000:]}')
    return str(csv)

completed = []
try:
    for job in JOBS:
        result_dir = OUT / 'resultados' / job['id']
        if (result_dir / 'completo.json').exists():
            completed.append(job['id'])
            print('SKIP COMPLETE', job['id'], flush=True)
            continue
        # Attempt labels also select the matching depth threshold in the summary script.
        attempts = (2 if job['final_settings']['sampling']==2000 else 1,) if args.final_settings else (1,2)
        for attempt in attempts:
            status(state='sampling', current=job['id'], attempt=attempt,
                   completed=completed, total=len(JOBS))
            print('START',job['id'],'attempt',attempt,flush=True)
            t0=time.monotonic()
            chain_files = [Path(job['dir']) / f'intento{attempt}_cadena{c}.csv' for c in range(1,5)]
            finished = all(p.exists() and 'Elapsed Time:' in p.read_text()[-1500:] for p in chain_files)
            if not finished:
                with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
                    list(pool.map(lambda c: chain(job,attempt,c), range(1,5)))
            print('SAMPLED',job['id'],'seconds',round(time.monotonic()-t0,1),flush=True)
            status(state='summarizing',current=job['id'],attempt=attempt,
                   completed=completed,total=len(JOBS))
            log = Path(job['dir']) / f'intento{attempt}_resumen.log'
            summary_script = Path(job['dir']) / f'intento{attempt}_resumir_ajuste.R'
            shutil.copy2(OUT/'codigo/resumir_ajuste.R', summary_script)
            with log.open('w') as stream:
                proc = subprocess.run(['Rscript','--vanilla',str(summary_script),
                                       job['id'],str(attempt)] + (['final'] if args.final_settings else []), stdout=stream,stderr=subprocess.STDOUT)
            if proc.returncode:
                raise RuntimeError(f'Summary failed: {log}\n{log.read_text()[-3000:]}')
            diagnostic=json.loads((Path(job['dir'])/'diagnostico_actual.json').read_text())
            if diagnostic['pass'] or attempt == 2 or args.final_settings:
                break
            print('RETRY',job['id'],diagnostic['diagnostics'],flush=True)
        completed.append(job['id'])
        print('DONE',job['id'],flush=True)
    status(state='complete',completed=completed,total=len(JOBS))
except Exception as error:
    status(state='error',current=job['id'],completed=completed,total=len(JOBS),error=str(error))
    raise
