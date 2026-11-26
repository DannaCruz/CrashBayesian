# Bayesian spectral analysis on road networks

This package contains models, aggregate inputs and computational results for road-segment crash counts in Barcelona and Bogotá. It evaluates the role of network topology, recorded road attributes and segment length across 20 fitted specifications.

## Reproduce the analysis

Start with the [finding-by-finding reproduction guide](REPRODUCIBILITY.md). It maps Tables 1–4, Figures 1–5, supplementary tables/figures and narrative metrics to the data, scripts, parameters and outputs. It distinguishes model replay from saved-result export and lists the unavailable source-event and historical simulation stages.

For a tabular check without sampling, use the [CSV-to-LaTeX table exporter](analysis/export_tables.py). The [validation record](analysis/reproduction_validation.json) states what was actually checked.

## Explore the results

- [Technical report](documents/manuscript.pdf)
- [Additional methods and diagnostics](documents/supplement.pdf)
- [Model comparison](analysis/reference/comparacion_modelos.csv)
- [Diagnostics for all 20 fits](analysis/reference/diagnosticos_20_especificaciones.csv)
- [Segment-length sensitivity](analysis/reference/sensibilidad_longitud_bayesiana.csv)

The separate [node–edge analysis](../Supplementary_Code/README.md) uses different inputs and model specifications.

## Model design

The main comparison combines four model families, two cities and two attribute specifications, giving 16 fits. Four additional NB/mixture fits estimate the length coefficient. The retained networks contain 2,540 segments in Barcelona and 6,079 in Bogotá.

Sixteen fits passed the joint sampler criterion. The four BYM2-type fits retained energy-exploration warnings after longer runs; some also retained ESS or R-hat warnings. Interpret numerical rankings together with these diagnostics.

The simulation material in the additional document uses an earlier model specification. It was not rerun for the ordered prior and does not validate the attribute-adjusted or estimated-length models.

## Package contents

- `documents/`: technical reports, editable LaTeX sources, bibliography, figures and tables. Compile here with `latexmk -pdf manuscript.tex` and `latexmk -pdf supplement.tex`.
- `analysis/inputs/`: aggregate counts, length offsets, design matrices, candidate spectral columns and binary segment adjacency.
- `analysis/models/`: Stan model source.
- `analysis/runs.json` and `analysis/initializations/`: final seeds, sampling settings and initializations.
- `analysis/reference/`: saved numerical summaries.
- `analysis/workflow/`: preparation, fitting, diagnostics and summary scripts. These scripts retain their `outputs/modelos_espectrales` and `work/modelos_espectrales` working layout.
- `analysis/run_model.py`: portable replay of a final fit from aggregate inputs.

## Check the inputs

From `analysis/`, run:

```sh
python3 run_model.py --check
```

This validates the aggregate inputs and settings without sampling.

## Run a model

Use an existing CmdStan installation and a C++ toolchain. The reference fits used CmdStan 2.36.0.

```sh
python3 run_model.py \
  --id Barcelona__Seleccion_espectral__con_covariables__longitud_libre \
  --cmdstan /path/to/cmdstan-2.36.0 \
  --output reruns/barcelona_length \
  --workers 4
```

Choose an ID from `runs.json` and a new output directory. The runner refuses an existing output directory. Seed and chain ID are passed separately to CmdStan. Compiler and platform differences can change MCMC trajectories.

The portable runner was input-checked; all 20 fits were not rerun during packaging. Completing a run does not establish convergence or predictive adequacy.

## Geographic preparation and figures

The optional preparation script requires `CRASHBAYESIAN_SOURCE_DIR` to identify the original network RData directory. Geographic preprocessing and map regeneration require the external source layers.

The figure/table generator accepts explicit source and destination directories. From this repository's root:

```sh
Rscript spectral-analysis/documents/generate_figures_tables.R \
  /path/to/modelos_espectrales spectral-analysis/documents
```

The source directory must contain the prepared geometry-bearing inputs and model summaries expected by the generator. It cannot regenerate maps from aggregate JSON inputs alone.

## Data scope and attribution

Individual crash records, geographic snapshots, original network RData files and complete original chain CSVs are not included. The shared model inputs start from prepared segment counts. Exact correspondence with dated municipal source records and temporal alignment of road attributes remain to be reconciled.

Road-network information derives from [OpenStreetMap contributors](https://www.openstreetmap.org/copyright). The associated source portals are [Open Data BCN](https://opendata-ajuntament.barcelona.cat/) and [Bogotá Open Data](https://datosabiertos.bogota.gov.co/). Existing source-data attribution and licensing remain applicable.

Outstanding checks include the Bogotá observation period, influential PSIS observations, residual dependence, zero-count calibration, sensitivity to the spectral basis and spatially blocked validation.

`MANIFEST_SHA256.json` records file hashes. Preparation-time code hashes in the reference results describe the original computation and may differ from later path-configured copies.
