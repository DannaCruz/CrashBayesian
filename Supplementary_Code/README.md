# Supplementary Code: Bogotá node–edge crash models

This archive starts from the exact aggregate analysis inputs used in the manuscript. It is separate from the earlier public network-base example. The latter uses a 2023 segment assignment and the principal connected component and is NOT the full 2007–2024 two-support analysis.

## Contents

- `runs.json`: 17 fit specifications, comprising four final primary CAR executions, six speed-sensitivity fits, five radius/length fits, and two original non-spatial Bayesian baselines. Four chains per fit; chain seed = base seed + chain id. Fixed adaptation constants and initialization files are included. Primary fits use 2,000 warm-up + 2,000 retained iterations; historical non-spatial baselines use 500 + 500. Baseline settings are documented by the archived shell workflow, not re-audited from original chains here.
- `inputs/`: response counts, standardized covariates, offsets, graph structure, fixed tuning constants and initializations in Stan JSON. Temporal inputs contain 13-year training and five-year test counts by unit, design matrices and geometric offsets. No individual crash records or coordinates are included.
- `models/`: original CAR, equivalent adaptive CAR with Jacobian, and non-spatial NB Stan models.
- `scripts/run_fit.py`: portable opt-in model runner, refusing existing output directories. `diagnose.R`, `diagnose_lp.R`, `compare_runs.R`: final-fit diagnostics and comparisons.
- `scripts/reproduce_simulation.R`: 250 replicates with original master seed 20260824 and sequential RNG consumption; there are no separately assigned per-replicate seeds. `original/04_support_simulation.R` preserves the source implementation. Its `safe_nb` helper falls back to quasi-Poisson on a fitting error; this behavior has not been silently removed.
- `scripts/reproduce_temporal.R`: the same non-spatial ML fits and 5/13 prediction scaling as the original, from aggregate inputs. `original/03_temporal_validation.R` documents the original event-to-period aggregation and requires the external processed geographic source objects.
- `reference/`: original saved simulation replicates, temporal metrics/predictions, baseline summaries, final CAR diagnostics/stability tables, coding tables and sensitivity results. Files under `reference/original` include superseded initial CAR summaries for provenance; use `reference/primary` for the revised final CAR results. No original-CAR output should replace the final fits.

## Quick verification (no Stan sampling)

From this directory, with Python 3 and R/MASS available:

```
python3 scripts/run_fit.py --check
Rscript scripts/reproduce_simulation.R NEW_SIMULATION_OUTPUT
Rscript scripts/reproduce_temporal.R NEW_TEMPORAL_OUTPUT
```

Both R commands refuse existing output directories and write into `results/` below the new directory. Compare with `reference/original/`. `REPRODUCTION_CHECK.csv` records validation in the supplied environment. All 250 saved simulation rows and summaries were reproduced exactly. Temporal coefficients agreed within 1e-15 and predictions within 5.5e-13. A raw Spearman difference of 1.03e-5 for intersection covariate predictions resulted from machine-precision changes to ties; rounding both saved and rerun predictions to 8, 10 or 12 decimals gives the same correlation. MAE, RMSE and capture are unchanged. The manuscript retains the archived metrics, not a replacement rerun.

## Optional long fits

Install CmdStan (the final primary chains used 2.36.0), a C++ toolchain, R packages cmdstanr and posterior, and set CMDSTAN to its installation directory. From this directory:

```
python3 scripts/run_fit.py --label segment_1 --output NEW_SEGMENT_RUN
```

Use a label from runs.json. The runner compiles locally, runs the four chains serially with the specified distinct seeds, and produces diagnostics for CAR fits. Completion does not imply convergence. Large fits need substantial memory, storage and runtime. The new runner was input-checked and syntax-checked; the 17 fits were not rerun for packaging. Machine/compiler differences can prevent bitwise-identical MCMC output.

## Reproduction boundary

This is an aggregate-input model-reproduction archive, not a complete raw-data-to-map distribution. Full posterior chains, geographic layer snapshots and individual crash data are not in this ZIP. Geographic preprocessing and re-rendering maps require the original external layers and assignment inputs. Exact model inputs and aggregate geometric audits are included so the model analysis does not depend on reconstructing those layers. The earlier public repository must not be substituted for them. Saved diagnostics cannot regenerate posterior traces without chains; the manuscript source package supplies the published figures. Historical and pilot inputs may be present, but only runs.json labels designate reported final fits.

Data attribution is in DATA_LICENSE.md; existing code licensing is in LICENSE. Publication of this prepared archive has not yet occurred; no public DOI or release is claimed.
