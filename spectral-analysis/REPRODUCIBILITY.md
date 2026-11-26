# Reproducing the spectral analysis

This guide maps each reported table, figure and numerical finding to its inputs, computation and output. Paths below are relative to `spectral-analysis/`, unless a command explicitly changes directory.

## Shared materials

- [Download all files currently shared in this repository (ZIP)](https://github.com/DannaCruz/CrashBayesian/archive/refs/heads/main.zip)
- [Aggregate analysis inputs and adjacency](analysis/inputs/)
- [Stan model](analysis/models/)
- [Portable model replay](analysis/run_model.py) and [final sampling configurations](analysis/runs.json)
- [Workflow scripts](analysis/workflow/) and [table exporter](analysis/export_tables.py)
- [Numerical reference outputs](analysis/reference/) and [validation record](analysis/reproduction_validation.json)
- [Technical documents](documents/)

The ZIP contains the files currently deposited; it does not add the missing original-event or historical spectral-simulation inputs listed below. All instructions and the finding-by-finding map follow on this page.

## 1. Reproduction scope

There are three different operations:

1. **Inspect saved results:** read the numerical files in `analysis/reference/`. This checks the reported summaries but does not refit models.
2. **Replay a model:** run `analysis/run_model.py` with the shared aggregate inputs, archived initialization, seed and final sampling settings. This produces new chains, rather than just plotting saved results.
3. **Reconstruct geographic inputs and downstream summaries:** use the workflow in `analysis/workflow/` with the external network RData inputs. This includes preprocessing, fitting, PSIS-LOO, posterior predictive checks and graph diagnostics. Map generation additionally requires the geometry-bearing prepared RDS objects.

**The complete original-event-to-finding workflow is not yet available.** Individual municipal event records, their exact source versions, the initial cleaning/assignment code, geographic source snapshots and complete original posterior chains are not in this repository. Network RData inputs are external. The supplied preparation script starts from those network files; it does not create them from original accident events. The matching historical spectral simulation code and per-replicate results are also unavailable here. These gaps are described below rather than represented as completed steps.

Full model reruns were not executed while preparing this guide. Saved diagnostics and an input check are not evidence of a successful fresh end-to-end run.

## 2. Software and fixed settings

For model replay, use Python 3, an installed CmdStan and a C++ toolchain. Reference fits used **CmdStan 2.36.0**. For the R workflow, install `sf`, `sfnetworks`, `igraph`, `Matrix`, `RSpectra`, `cmdstanr`, `posterior`, `loo`, `statmod`, `matrixStats`, `jsonlite`, `dplyr` and `ggplot2`. See `analysis/reference/sesion_R.txt` and the package-version information supplied with the computational archive. LaTeX compilation requires pdfLaTeX, BibTeX and the packages in `documents/preamble.tex`.

- 16 fits: two cities × four model families × two attribute settings.
- Four additional NB/mixture fits estimate length, giving 20 total.
- `model_type`: 0 NB; 1 dense spectral; 2 continuous-mixture spectral; 3 spectral BYM2.
- U = no attributes, unit length coefficient; C = attributes, unit length coefficient; L = attributes, estimated length coefficient.
- Each fit uses four chains, with chain IDs 1–4 and the base seed in `analysis/runs.json`. The base seed and chain ID are passed separately to CmdStan.
- Final BYM2 settings: 2,000 warm-up + 2,000 retained iterations per chain, `adapt_delta=0.995`, maximum depth 14.
- Other final fits: 1,000 + 1,000, `adapt_delta=0.99`, maximum depth 13.
- Initializations are in `analysis/initializations/`. Do not replace them if replaying archived final settings.

All shared data paths in the portable runner are resolved from its own location. Output directories and the CmdStan installation are command-line arguments; no source editing is required. The older R workflow intentionally retains an `outputs/modelos_espectrales` and `work/modelos_espectrales` working layout; section 5 constructs it explicitly.

The delivered analysis is scripted. No GUI procedure is required for these commands, so GUI screenshots are not applicable to this workflow. If any missing upstream stage used GUI software, its procedure, screenshots and exact parameters still need to accompany that stage.

## 3. Check and replay the shared aggregate inputs

From `spectral-analysis/`:

```sh
python3 analysis/run_model.py --check
```

This validates dimensions, finite values, nonnegative counts and initialization compatibility for all 20 specifications. It does not sample or reconstruct original events.

Replay one final fit:

```sh
python3 analysis/run_model.py \
  --id Barcelona__Seleccion_espectral__con_covariables__longitud_libre \
  --cmdstan /path/to/cmdstan-2.36.0 \
  --output reruns/barcelona_length \
  --workers 4
```

`--cmdstan` identifies a software installation, not a data path. Choose a new output directory; the runner refuses an existing one. It writes four `chain_*.csv` files, their logs, `data.json`, `initialization.json` and `run.json`. Compiler/platform differences can change MCMC trajectories. Completing sampling does not establish convergence.

To replay all final configurations without editing code:

```sh
python3 - <<'PY'
from pathlib import Path
import json, os, subprocess, sys
root = Path('analysis').resolve()
cmdstan = os.environ.get('CMDSTAN')
if not cmdstan:
    raise SystemExit('Set CMDSTAN to the installed CmdStan directory first.')
output = Path('reruns/all_final_fits').resolve()
if output.exists():
    raise SystemExit('Choose a fresh checkout/output location; output already exists.')
runs = json.loads((root / 'runs.json').read_text())
for run in runs:
    subprocess.run([
        sys.executable, str(root / 'run_model.py'), '--id', run['id'],
        '--cmdstan', cmdstan, '--output', str(output / run['id']),
        '--workers', '4'
    ], check=True)
PY
```

Set `CMDSTAN` to your installation before that command. New chains must then be summarized and assessed; they are not table/figure outputs by themselves. The archived summary implementation and its expected layout are described in sections 5–6. The portable `run_model.py` does not automatically adapt its chain filenames/layout to that older summary workflow.

## 4. Original records to network inputs: missing stage

For each city, a complete reconstruction requires:

1. Exact municipal dataset version or immutable snapshot, download query, date range, license and checksum.
2. Event inclusion/exclusion criteria, duplicate handling, coordinate checks and CRS transformations.
3. Dated OSM/network extraction, road boundary, segmentation, parallel-road handling and geometry cleaning.
4. Event-to-segment assignment code: tolerances, intersection cases, unmatched events and count aggregation.
5. The transformation that produced `roads$Y = log(1 + count)` and the code saving the network RData objects.
6. Reconciliation of the Bogotá observation period and the incompatible `Y_edges` object. Do not infer this period from the separate node–edge analysis elsewhere in the repository.

Expected network files are `Eixample_network_ready_for_models.RData` and `Bogota_network_ready_for_models.RData`. The portable aggregate JSON files record the corresponding source basenames and MD5 hashes. A current portal download is not an interchangeable replacement for an unidentified historical snapshot.

`analysis/workflow/preparar_datos.R` loads those network files and starts **after** the missing steps above. Recovering integer counts with `round(expm1(roads$Y))` verifies a numerical transformation, not deduplication or original event assignment.

## 5. Geographic workflow with external network inputs

This route requires the two external network RData files. It cannot be run from aggregate JSON alone. From `spectral-analysis/`, create a fresh working layout:

```sh
mkdir -p workspace/outputs/modelos_espectrales/codigo
mkdir -p workspace/outputs/modelos_espectrales/datos
mkdir -p workspace/outputs/modelos_espectrales/resultados
mkdir -p workspace/outputs/modelos_espectrales/figuras
mkdir -p workspace/work/modelos_espectrales
mkdir -p workspace/data/originals
cp analysis/workflow/* workspace/outputs/modelos_espectrales/codigo/
```

Place the two network RData files in `workspace/data/originals/`. Do not create substitute source records. Then:

```sh
cd workspace
export CRASHBAYESIAN_SOURCE_DIR=data/originals
Rscript --vanilla outputs/modelos_espectrales/codigo/preparar_datos.R
Rscript --vanilla outputs/modelos_espectrales/codigo/crear_trabajos.R
Rscript --vanilla outputs/modelos_espectrales/codigo/agregar_sensibilidad_longitud.R
Rscript --vanilla outputs/modelos_espectrales/codigo/diagnosticar_longitud.R
python3 outputs/modelos_espectrales/codigo/compilar_modelo.py --cmdstan /path/to/cmdstan-2.36.0
python3 outputs/modelos_espectrales/codigo/correr_modelos.py --workers 4
Rscript --vanilla outputs/modelos_espectrales/codigo/compilar_resultados.R
Rscript --vanilla outputs/modelos_espectrales/codigo/cerrar_corrida.R
```

This older executor uses generic starts and its original diagnostic-triggered retry rule, so it is a **statistical rerun**, not a replay of every archived final initialization. For archived final settings, use section 3. The Bogotá BYM2 initialization helper documents a specialized historical start; do not treat a different start as bitwise-identical historical replay.

Preprocessing uses seed `20261120`, checks geometry/graph order, preserves segment order and parallel features, selects the component with the greatest number of segments, constructs binary shared-endpoint adjacency, and computes 61 leading normalized-adjacency eigenpairs. Signs are fixed deterministically. Perron is omitted; 20 columns are centered and RMS-scaled. Length is in kilometres. Attribute references are residential roads, two lanes and asphalt; invalid lanes are coded unknown. Outputs include:

- `datos/*_preparado.rds` and `datos/*_segmentos.csv`;
- `resultados/auditoria_datos.csv`, `covariables_niveles.csv`, and `*_espectro.csv`;
- per-fit chain inputs, chains/logs and numerical results.

An independent local check reran this preprocessing from the recorded network RData files and matched the archived audit counts, category counts, lengths, eigenvalues and eigen-equation residuals. This did not reconstruct municipal source events.

The scripts are commented at the alignment checks, eigenvector sign convention, mixture summaries, integrated likelihood, quadrature and predictive checks. Read their comments alongside the stage descriptions below.

## 6. Posterior, predictive and display computation

`resumir_ajuste.R` is the numerical implementation. It reads the four chains in the older workflow layout and:

1. Computes rank-normalized R-hat, bulk/tail ESS, divergences, depth hits and chain E-BFMI; BYM2 also checks independent effects.
2. Computes coefficient summaries and exponentiated count ratios, mixture membership and the posterior Poisson-binomial distribution of K.
3. Reconstructs signed contributions and the complete structured field using all coefficients.
4. Builds likelihood matrices from up to 500 evenly spaced draws per chain, then computes PSIS-LOO and its diagnostics.
5. For BYM2 LOO, integrates the held-out independent effect with 30-node Gaussian quadrature, checks 60 nodes and switches to 60 if the maximum check difference exceeds 0.001.
6. Uses 400 predictive replicates with seed `job_seed + 40000`; computes zero fraction, variance, maximum, 99th percentile and mean, and 100 selected residual Moran comparisons.

The original joint sampler criterion is no divergences/depth hits, maximum R-hat ≤1.01, minimum bulk ESS ≥400 and minimum chain E-BFMI ≥0.3, plus applicable BYM2 local-effect checks. A flagged fit is not established as converged.

For geographic plot/table regeneration, return to `spectral-analysis/` after section 5 and run:

```sh
Rscript --vanilla documents/generate_figures_tables.R \
  workspace/outputs/modelos_espectrales documents
```

This needs geometry-bearing RDS inputs and the numerical outputs. Aggregate JSON does not contain map geometries. Saving plots from reference outputs alone is display regeneration, not model reproduction. A separate `analysis/export_tables.py` command below exports tabular findings directly from saved CSV summaries and makes this distinction explicit.

Compile the documents from `documents/`:

```sh
latexmk -pdf manuscript.tex
latexmk -pdf supplement.tex
```

Narrative estimates are not automatically rewritten after a statistical rerun. Compare every mapped numerical finding with the new outputs before interpreting them.

## 7. Finding-by-finding instructions

The common computation for model-based findings is sections 3 or 5, followed by section 6. The paths in this section identify the exact outputs and required transformations, rather than only a final plotted image.

Fit IDs have the form `CITY__MODEL__SPECIFICATION`:

- CITY: `Barcelona` or `Bogota`.
- MODEL: `NB`, `CAR_espectral`, `Seleccion_espectral`, `BYM2_espectral`.
- SPECIFICATION: `sin_covariables` (U), `con_covariables` (C), or `con_covariables__longitud_libre` (L).
- L exists only for NB and mixture models.

Saved per-fit files are in `analysis/reference/fits/FIT_ID/`. Recomputed workflow files are in `workspace/outputs/modelos_espectrales/resultados/FIT_ID/`.

### Main tables

| Finding | Required upstream computation and output |
|---|---|
| Table 1: retained analytical samples | Original events → external network inputs → `preparar_datos.R` → `auditoria_datos.csv`. Read retained segment/intersection counts, recovered total, represented kilometres and zero fraction. Source-event reconciliation remains missing. |
| Table 2: posterior mixture membership | Fit the six mixture U/C/L specifications → `resumir_ajuste.R` → `seleccion_modos.csv` and `distribucion_K.csv`. Sum PIP for E(K), take cumulative 5th/95th K quantiles and list modes with PIP >0.5. |
| Table 3: conditional count ratios | Fit the two mixture L models → `covariables.csv`. Use exponentiated median, 5th and 95th coefficient percentiles for the six substantive contrasts; retain unknown-category coefficients in estimation. |
| Table 4: all 20 PSIS-LOO results | All 20 fits → integrated/ordinary likelihood calculations in `resumir_ajuste.R` → `prediccion_loo.csv` and `diagnosticos.csv`. Report ELPD, ELPD SE, Pareto exceedances and the joint sampling flag. |

To export these four table bodies from the saved CSVs, without sampling:

```sh
python3 analysis/export_tables.py --output derived/reference_tables
```

That exporter also writes S1–S5. It refuses an existing output directory. For newly computed workflow results:

```sh
python3 analysis/export_tables.py \
  --fits workspace/outputs/modelos_espectrales/resultados \
  --output derived/rerun_tables
```

This command starts at numerical summaries; the fitting and likelihood stages above are still required for a model rerun.

### Main figures

| Finding | Required computation and plotting input |
|---|---|
| Figure 1: leading eigenvalues | Source network → adjacency → normalized operator → leading 61 eigenpairs in `preparar_datos.R` → `Barcelona_espectro.csv`/`Bogota_espectro.csv` → `spectrum.pdf`. Only the upper spectrum is computed. |
| Figure 2: first three raw candidate eigenvectors | Same preprocessing → prepared `U_raw` and segment geometry. Map columns 1–3 multiplied by sqrt(n), without centering → `raw_network_modes.pdf`. Geometries and raw eigenvectors are required. |
| Figure 3: U/C/L membership | Six mixture fits → `seleccion_modos.csv` → `compilar_resultados.R` → `seleccion_con_sensibilidad_longitud.csv` → `mode_membership.pdf`. Keep the same candidate basis within each city. |
| Figure 4: signed contributions | Two mixture L fits → `contribuciones_modos.csv`. Rank RMS contributions; multiply each selected mean coefficient by its matching B column and map the three largest with retained geometry → `adjusted_contributions.pdf`. |
| Figure 5: complete fields | Two mixture L fits → `resumen_segmentos.csv`, `structured_field_mean`. Map the field combining all 20 coefficients, with geometry → `adjusted_field.pdf`. Do not threshold the fitted field at PIP >0.5. |

For Figures 1–5, run the geometric display command in section 6 after the corresponding upstream calculations. Figures 2, 4 and 5 cannot be recreated as geographic maps with the aggregate JSON alone.

### Empirical supplementary findings

The compact supplementary document contains S1–S5 and Figure S1; these empirical labels also occur at the start of the older extended supplementary document stored in this repository.

| Finding | Required upstream computation and output |
|---|---|
| Table S1: category counts | Preprocessing → `covariables_niveles.csv`; count highway, lane and surface factor levels on retained segments. |
| Table S2: all count ratios | Two mixture L fits → `covariables.csv`; include unknown-category coefficients and their marginal 90% intervals. |
| Table S3: sampler diagnostics | Four chains for each of 20 fits → `diagnosticos.csv`; draws/chain, maximum R-hat, minimum bulk ESS and minimum chain E-BFMI, with BYM2 local-effect checks. |
| Table S4: predictive checks | Two mixture L fits → 400 negative-binomial replications in `resumir_ajuste.R` → `comprobaciones_predictivas.csv`; observed zero fraction, variance, maximum, q99 and mean, replicated median/90% interval and upper-tail probability. |
| Table S5: residual Moran I | NB and mixture L fits for both cities → centered Pearson residuals and binary adjacency → `autocorrelacion_residual.csv`; observed median, replicated median/interval and upper-tail probability. |
| Figure S1: field interval sign | Two mixture L fits → `resumen_segmentos.csv`; classify q05 >0, q95 <0 or interval containing zero, join geometry → `field_uncertainty.pdf`. Intervals are marginal, not a simultaneous hotspot test. |

### Metrics reported in the text

- **E(K), individual PIP and thresholded mode lists:** six mixture `seleccion_modos.csv` files. This includes Barcelona 6.92/2.27/1.31 and Bogotá 4.44/2.89/2.47, highest PIP and modes above 0.5. K intervals require `distribucion_K.csv`, not quantiles of the draw-wise probability sum.
- **Length elasticities and 90% intervals:** four L `elasticidad_longitud.csv` files. The total elasticity is 1 plus the `log_length_deviation` coefficient. The approximate 0.59 estimates are conditional associations, not traffic-volume-adjusted effects.
- **ELPD gains and paired SEs:** use matched `loo.rds` pointwise `elpd_loo` vectors on identical segments. Gain = sum of differences; paired SE = sqrt(n × variance of those differences). `compilar_resultados.R` computes the length and attribute comparisons. Mixture-vs-NB L contrasts use those same paired vectors. Do not subtract aggregate ELPD SEs. Pointwise LOO objects are regenerated from chains; they are not included among the public reference CSVs.
- **Reported gains:** mixture length gains 76.3/265.3, NB length gains 76.8/252.3, mixture attribute gains 36.6/404.7 and mixture-vs-NB L gains 18.1/95.4 follow the preceding paired calculation. Saved comparisons are in `ganancia_covariables.csv` and `sensibilidad_longitud_bayesiana.csv`.
- **Pareto k, including maximum 0.743:** per-fit `prediccion_loo.csv` and regenerated `loo.rds`. The threshold is min(1 − 1/log10(S), 0.7), with S the likelihood draw count.
- **R-hat, ESS, E-BFMI, divergences and depth hits:** per-fit `diagnosticos.csv`. The reported BYM2 E-BFMI range 0.191–0.277 and ESS/R-hat warnings must remain attached to their score comparisons.
- **Predictive means, variances, zeros, maxima and q99:** `comprobaciones_predictivas.csv` for the stated city/model/specification. The numerical statements in the results use the same workflow as S4; C-model comparisons use C files.
- **Residual Moran I under C and L:** the corresponding `autocorrelacion_residual.csv` files, with centered Pearson residuals and the same binary adjacency.
- **Original/excluded counts, recovered totals, missing categories and invalid lane=22:** `auditoria_datos.csv`, `covariables_niveles.csv` and original network attribute fields. Original event deduplication/assignment is still outside this workflow.
- **Eigen-equation residual below 6e-11:** residual column in both `*_espectro.csv` files and `eigen_residual` in prepared RDS objects. Geometry/graph order checks are in preprocessing. A separate independent topology-verification script was not supplied.

### Historical supplementary material

The older extended supplement additionally includes the following labels. Their generating code and per-replicate outputs were not located, so reproducing the reported historical values is **blocked**:

| Historical finding | Inputs/code required but not shared |
|---|---|
| Table S6: composite difficulty settings | Exact spectral experiment configuration and generator. The table describes settings but does not execute them. |
| Table S7: pooled simulation performance | Outputs from 1,200 datasets and four Bayesian fits each; per-replicate RMSE, predictive scores, coverage and paired differences. |
| Table S8: mode recovery | Per-replicate generating active sets and posterior mode summaries for sparse truth; AUC, TPR, FDR and E(K). |
| Table S9: simulation computation | All 4,800 fit logs, elapsed times, divergences, depth/ESS summaries and exact pass-classification implementation. |
| Figure S2: RMSE trajectories | Per-cell/per-replicate RMSE outputs and plotting code (`Fig04.png` is only the saved image). |
| Figure S3: relative predictive scores | Per-replicate test-score differences and Monte Carlo SE calculations (`Fig05.png`). |
| Figure S4: mode-recovery trajectories | Per-replicate AUC/TPR/FDR and aggregation/plotting code (`Fig06.png`). |
| Figure S5: Barcelona representative fields | Matching historical network/basis, generator, model-specific seeds, moderate difficulty replicate 25 and plotting code (`Fig16.png`). |
| Figure S6: Bogotá representative fields | Corresponding Bogotá inputs and generator for replicate 25 (`Fig17.png`). |

The described design is 2 networks × 4 mechanisms × 3 difficulty settings × 50 replications = 1,200 datasets and 4,800 Bayesian fits. The documented task seed is `910000 + 1000*scenario_id + replicate`; the exact scenario ordering and model-specific seed derivation must be supplied, not guessed. Historical performance claims, coverage, paired Monte Carlo SEs, pass rates and execution counts cannot be reconstructed from aggregate tables/images alone.

`Supplementary_Code/scripts/reproduce_simulation.R` belongs to a **different**, 250-replicate node–edge-support experiment. It must not be substituted for this spectral experiment.

## 8. Validation record and remaining work

Verified during preparation of the accompanying computational package:

- Local network preprocessing matched archived numerical audits and spectra.
- Final settings and 20 initializations matched their recorded SHA-256 hashes.
- The working-package figures and nine empirical table bodies were regenerated.
- The local compact LaTeX supplement compiled to three pages; the extended repository PDF is a different document version.

The public table exporter was checked against the supplied per-fit/reference CSVs. These checks do not constitute a fresh MCMC run. A complete reproducibility record still requires original-event provenance/processing, matching historical spectral simulation code if those findings are retained, full chains/logs and a successful clean end-to-end execution with every numerical comparison checked. Keep scientific limitations and diagnostic flags when reporting reruns.
