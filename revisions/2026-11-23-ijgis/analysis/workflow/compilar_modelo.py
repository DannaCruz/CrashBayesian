"""Compile the delivered Stan source locally with an existing CmdStan installation."""
from pathlib import Path
import argparse
import shutil
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument('--cmdstan', type=Path,
                    default=Path.home()/'.cmdstan/cmdstan-2.36.0')
args = parser.parse_args()
source = Path('outputs/modelos_espectrales/codigo/modelo_nb_espectral.stan').resolve()
target_dir = Path('work/modelos_espectrales/stan').resolve()
target_dir.mkdir(parents=True, exist_ok=True)
target = target_dir / source.name
shutil.copy2(source, target)
subprocess.run(['make', '-C', str(args.cmdstan), str(target.with_suffix('')),
                'PRECOMPILED_HEADERS=false'], check=True)
