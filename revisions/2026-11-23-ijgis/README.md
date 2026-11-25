# IJGIS manuscript revision: November 2026

This directory records the revised manuscript **Bayesian spectral selection for road-segment crash counts: road attributes and length sensitivity** and its 20 fitted empirical specifications. It is separate from the repository's existing [17-specification supplementary archive](../../Supplementary_Code/README.md).

## Read the documents

- [Main manuscript (PDF)](manuscript/manuscript.pdf)
- [Supplementary material (PDF)](manuscript/supplement.pdf)
- [Editable manuscript](manuscript/manuscript.tex), [supplement](manuscript/supplement.tex), and [bibliography](manuscript/references.bib)
- [Model comparison](analysis/reference/comparacion_modelos.csv)
- [All 20 sampler diagnostics](analysis/reference/diagnosticos_20_especificaciones.csv)
- [Estimated-length sensitivity](analysis/reference/sensibilidad_longitud_bayesiana.csv)

The supplement contains the complete simulation material retained from the preceding manuscript, including its design, tables and figures. Those simulations were **not rerun** for the revised ordered prior and do not constitute validation of the revised specification. No new simulations are part of this publication update.

## Empirical scope

The principal comparison contains four model families (NB, dense spectral, continuous-mixture spectral and spectral BYM2), two cities and two attribute specifications: 16 fits. Four additional NB/mixture fits estimate the length coefficient, giving 20 specifications in total. The retained supports have 2,540 segments in Barcelona and 6,079 in Bogotá. Barcelona uses the `Eixample_network_ready_for_models.RData` source, not the 447-segment subnetwork.

Sixteen fits passed the stated joint sampler criterion. All four BYM2 fits retained energy-exploration warnings after longer runs; some also retained ESS or R-hat warnings. Scores and saved output must be interpreted with these diagnostics. The main manuscript is an **author-review revision**, not a declaration that all scientific checks for submission have been completed. Event provenance, the Bogotá observation period, temporal alignment of attributes, remaining predictive discrepancies, influential LOO cases and blocked spatial validation remain documented limitations.

## Files and reproduction boundary

- `manuscript/`: both PDFs, LaTeX sources, bibliography, English figures, retained simulation figures and source-linked tables. Compile inside this directory with `latexmk -pdf manuscript.tex` and `latexmk -pdf supplement.tex`.
- `analysis/inputs/`: aggregate counts, length offsets, design matrices, candidate spectral columns and binary segment adjacency. No individual crash records, coordinates, geographic layer snapshots or original RData files are included.
- `analysis/models/`: the Stan model used for the revised fits.
- `analysis/runs.json` and `analysis/initializations/`: final seeds, settings and initializations. Seed and chain id are passed separately to CmdStan, matching the recorded workflow.
- `analysis/reference/`: saved numerical summaries for the revised fits. Passing a sampler criterion does not imply distributional adequacy.
- `analysis/workflow/`: the analysis workflow used in the working project. Its relative `outputs/modelos_espectrales` and `work/modelos_espectrales` layout is retained. The optional preparation script accepts `CRASHBAYESIAN_SOURCE_DIR` instead of a personal filesystem path. Geographic preparation and map regeneration require the external source layers.
- `analysis/run_model.py`: portable, opt-in replay of one final fit from the aggregate inputs. This runner was input-checked, not used to rerun all fits during packaging. Machine/compiler differences can change MCMC trajectories.

From `analysis/`, verify the aggregate inputs without sampling:

```sh
python3 run_model.py --check
```

An optional model replay requires an existing CmdStan installation and C++ toolchain:

```sh
python3 run_model.py --id Barcelona__Seleccion_espectral__con_covariables__longitud_libre \
  --cmdstan /path/to/cmdstan-2.36.0 --output reruns/barcelona_length --workers 4
```

The runner refuses an existing output directory and does not run simulations. Complete original chain CSV files and geographic objects remain outside this public package. Derived field summaries use anonymous analysis segment indices, not identifiable event records.

## Attribution and evidence status

Road-network information derives from [OpenStreetMap contributors](https://www.openstreetmap.org/copyright); the underlying ODbL terms remain applicable. The manuscript identifies Open Data BCN and Bogotá's municipal open-data portal as source contexts for the prepared count datasets. Exact event-level correspondence and the relevant source versions remain to be reconciled; these aggregate inputs should not be represented as independently verified municipal event records. This directory makes no new blanket license grant over source data or manuscript content. Existing repository licensing and source attribution are preserved.

`MANIFEST_SHA256.json` records hashes for the files in this revision. Preparation-time code hashes in the reference results describe the original working run; they need not match the public copy of the path-configured preparation script or the new portable runner.
