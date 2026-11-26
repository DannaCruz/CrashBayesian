"""Replay one archived final fit; --check never samples or runs simulations."""
from pathlib import Path
import argparse
import concurrent.futures
import json
import math
import os
import shutil
import statistics
import subprocess

HERE = Path(__file__).resolve().parent

def model_data(run):
    city = json.loads((HERE / 'inputs' / f"{run['city']}.json").read_text())
    n = city['N']
    x = [list(row) for row in city['X']] if run['with_covariates'] else [[] for _ in range(n)]
    if run['free_length']:
        center = statistics.fmean(city['log_exposure'])
        x = [row + [offset - center] for row, offset in zip(x, city['log_exposure'])]
    spatial = run['model_type'] != 0
    return dict(N=n, P=len(x[0]), X=x, y=city['y'], log_exposure=city['log_exposure'],
                M=20 if spatial else 0, B=city['B'] if spatial else [[] for _ in range(n)],
                **{'lambda': city['lambda'] if spatial else []},
                model_type=run['model_type'], mixture=int(run['model_type'] == 2),
                iid=int(run['model_type'] == 3), gamma=city['gamma'],
                spike_sd=city['spike_sd'],
                field_normalizer=city['field_normalizer'] if spatial else 1)

def validate(run):
    d = model_data(run)
    assert len(d['y']) == len(d['X']) == len(d['B']) == len(d['log_exposure']) == d['N']
    assert all(isinstance(y, int) and y >= 0 for y in d['y'])
    assert all(len(x) == d['P'] and all(map(math.isfinite, x)) for x in d['X'])
    assert all(len(b) == d['M'] and all(map(math.isfinite, b)) for b in d['B'])
    assert len(d['lambda']) == d['M'] and all(-1 <= v <= 1 for v in d['lambda'])
    assert all(map(math.isfinite, d['log_exposure']))
    assert d['field_normalizer'] > 0
    init = json.loads((HERE / run['initialization']).read_text())
    assert len(init['beta']) == d['P'] and len(init['b']) == d['M']
    assert len(init['v_raw']) == (d['N'] if d['iid'] else 0)
    assert init['phi'] > 0 and run['chains'] == 4
    return d

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    parser.add_argument('--id', help='Exact id from runs.json')
    parser.add_argument('--cmdstan', type=Path, default=os.environ.get('CMDSTAN'))
    parser.add_argument('--output', type=Path)
    parser.add_argument('--workers', type=int, default=1)
    args = parser.parse_args()
    runs = json.loads((HERE / 'runs.json').read_text())
    for run in runs:
        validate(run)
    if args.check:
        print(f'Validated {len(runs)} final specifications. No sampling or simulations executed.')
        return
    if not args.id or not args.output or not args.cmdstan:
        parser.error('Use --check, or provide --id, --cmdstan and a new --output directory.')
    selected = [r for r in runs if r['id'] == args.id]
    if len(selected) != 1:
        parser.error('Unknown or ambiguous fit id.')
    if args.output.exists():
        parser.error('Output directory already exists; choose a new directory.')
    if not 1 <= args.workers <= 4:
        parser.error('--workers must be between 1 and 4.')
    run = selected[0]
    output = args.output.resolve()
    output.mkdir(parents=True)
    source = HERE / 'models' / 'modelo_nb_espectral.stan'
    build = output / 'stan'
    build.mkdir()
    target = build / source.name
    shutil.copy2(source, target)
    executable = target.with_suffix('')
    subprocess.run(['make', '-C', str(Path(args.cmdstan).resolve()), str(executable),
                    'PRECOMPILED_HEADERS=false'], check=True)
    data_file = output / 'data.json'
    data_file.write_text(json.dumps(model_data(run), allow_nan=False))
    init_file = output / 'initialization.json'
    shutil.copy2(HERE / run['initialization'], init_file)
    (output / 'run.json').write_text(json.dumps(run, indent=2))
    def chain(chain_id):
        command = [str(executable), 'sample', f"num_samples={run['sampling']}",
                   f"num_warmup={run['warmup']}", 'save_warmup=0', 'thin=1',
                   'adapt', f"delta={run['adapt_delta']}", 'algorithm=hmc', 'engine=nuts',
                   f"max_depth={run['max_depth']}", 'data', f'file={data_file}',
                   f'init={init_file}', 'random', f"seed={run['seed']}", f'id={chain_id}',
                   'output', f"file={output / f'chain_{chain_id}.csv'}", 'refresh=200', 'sig_figs=8']
        with (output / f'chain_{chain_id}.log').open('w') as stream:
            subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT, check=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        list(pool.map(chain, range(1, 5)))
    print(f'Completed {run["id"]}. Review sampler and predictive diagnostics before interpretation.')

if __name__ == '__main__':
    main()
