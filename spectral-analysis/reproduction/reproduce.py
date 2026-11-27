"""Reproduce the record-count analysis from shared source records and fixed roads.

Run from any directory. All project paths are resolved relative to this file.
Outputs must be new: supplied reference results are never used as fitted results.
"""
from pathlib import Path
import argparse, hashlib, json, shutil, subprocess, sys
from configurar_corrida_final import configure
ROOT=Path(__file__).resolve().parent
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--output',type=Path,required=True,help='New output directory')
p.add_argument('--cmdstan',type=Path,required=True,help='Installed CmdStan directory, reference version 2.36.0')
p.add_argument('--workers',type=int,default=4)
p.add_argument('--skip-pdf',action='store_true',help='Produce tables and figures without compiling LaTeX')
a=p.parse_args()
def run(cmd,cwd=ROOT):
 print('+',' '.join(map(str,cmd)),flush=True)
 subprocess.run(list(map(str,cmd)),cwd=cwd,check=True)
if a.output.exists():p.error('Output already exists; choose a new directory')
if a.workers<1:p.error('Workers must be positive')
cmdstan=a.cmdstan.resolve()
if not (cmdstan/'makefile').is_file():p.error('CmdStan makefile not found')
# Check all immutable source snapshots before processing any record.
for name,expected in json.loads((ROOT/'reproduction/input_hashes.json').read_text()).items():
 if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=expected:
  raise SystemExit('Source checksum mismatch: '+name)
run(['Rscript','--vanilla','-e',"pkgs<-c('sf','sfnetworks','igraph','Matrix','RSpectra','jsonlite','cmdstanr','posterior','loo','ggplot2','dplyr','tidyr','statmod');stopifnot(all(vapply(pkgs,requireNamespace,logical(1),quietly=TRUE)));cat(R.version.string,'\\n');print(vapply(pkgs,function(p)as.character(packageVersion(p)),character(1)))"])
out=a.output.resolve();out.mkdir(parents=True)
model=out/'modelos_espectrales'
for folder in ['codigo','inicializaciones_finales']:
 shutil.copytree(ROOT/'modelos_espectrales'/folder,model/folder)
shutil.copy2(ROOT/'modelos_espectrales/configuraciones_finales.json',model/'configuraciones_finales.json')
for folder in ['datos','resultados','figuras']:(model/folder).mkdir()
(out/'work/modelos_espectrales').mkdir(parents=True)
# Preserve original record units, repeats, period and original nearest-road rules.
run(['Rscript','--vanilla',ROOT/'reproduction/rebuild_from_sources.R','historical_records',ROOT/'source_data',out/'source_reconstruction'])
shutil.copytree(out/'source_reconstruction/networks',model/'datos/originales')
def r(script):run(['Rscript','--vanilla','modelos_espectrales/codigo/'+script],out)
r('preparar_datos.R');r('crear_trabajos.R');r('agregar_sensibilidad_longitud.R')
configure(model,out/'work/modelos_espectrales/trabajos.json')
r('diagnosticar_longitud.R')
run([sys.executable,'modelos_espectrales/codigo/compilar_modelo.py','--cmdstan',cmdstan],out)
# Four chains for each of 20 fits; seeds and budgets come from recorded settings.
run([sys.executable,'modelos_espectrales/codigo/correr_modelos.py','--final-settings','--workers',a.workers],out)
r('compilar_resultados.R');r('cerrar_corrida.R')
run(['Rscript','--vanilla',ROOT/'reproduction/collect_metrics.R',model,out/'metrics'])
shutil.copytree(ROOT/'documents',out/'documents',ignore=shutil.ignore_patterns('*.pdf','*.aux','*.log','*.bbl','*.blg','figures','tables'))
run([sys.executable,ROOT/'reproduction/export_tables.py','--fits',model/'resultados','--output',out/'documents/tables'])
run(['Rscript','--vanilla',ROOT/'documents/generate_figures_tables.R',model,out/'documents'])
if not a.skip_pdf:
 for stem in ['manuscript','manuscript_with_authors','supplement']:
  if (out/'documents'/(stem+'.tex')).is_file():
   run(['latexmk','-pdf','-interaction=nonstopmode','-halt-on-error',stem+'.tex'],out/'documents')
run([sys.executable,ROOT/'reproduction/verify_results.py','--fits',model/'resultados','--reference',ROOT/'modelos_espectrales/resultados','--output',out/'verification.json'])
print('Completed record-source-to-findings workflow:',out)
