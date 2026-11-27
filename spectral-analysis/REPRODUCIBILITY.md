# Reproduce the network record-count analysis

Download the [complete repository ZIP](https://github.com/DannaCruz/CrashBayesian/archive/refs/heads/main.zip), extract it and open `spectral-analysis/reproduction`. This directory includes the source-record projections, frozen road geometries, preprocessing, model, settings, initializations, reference summaries and document sources. No municipal portal login, external network RData, hand-edited path or graphical software procedure is required. Internet access is required only to install the software.

## Response units and immutable inputs

Barcelona uses accident-report **rows**, including repeated report identifiers, from the frozen 2023 Guardia Urbana file: 7,724 rows and 7,721 incident identifiers before geographical filtering. Bogotá uses **person-involvement rows**, including multiple people per crash, from the frozen 2007–2025 extract: 468,551 rows and 324,962 incident identifiers before filtering. These responses are not distinct-crash counts, comparable city incidence or traffic-adjusted risks.

`source_data/*_records.csv.gz` retains every original row and its order, incident identifier, year, source date and coordinates, excluding unused personal attributes. `source_projection_manifest.json` records original and projection hashes and row counts. `reproduction/project_source_columns.py` documents and implements this projection if the full original files are available; they are not needed to run the shared projection. No precomputed segment response is in these event-level files. The frozen road snapshots contain geometry and explanatory attributes, **no outcome**. `reproduction/input_hashes.json` verifies all shared source snapshots before execution.

Source attribution: Barcelona City Council, Guardia Urbana accident dataset, original `2023_accidents_gu_bcn.csv`, [dataset catalogue](https://datos.gob.es/es/catalogo/l01080193-accidentes-gestionados-por-la-guardia-urbana-en-la-ciudad-de-barcelona), CC BY 4.0. Bogotá Secretaría Distrital de Movilidad, original `accidentes_bogota.csv.gz`, public SIMUR accident service, [source FeatureServer](https://sig.simur.gov.co/arcgis/rest/services/Accidentalidad/WSAcidentalidad_Publico/FeatureServer/2), [municipal catalogue](https://datosabiertos.bogota.gov.co/dataset/historico-siniestros-bogota-d-c), CC BY 4.0. The archived extract, rather than a live query, defines the observations. Road geometry and attributes: © OpenStreetMap contributors, [ODbL](https://www.openstreetmap.org/copyright). Exact original OSM query dates were not recovered; frozen snapshots define the model support and preserve geometry and road-feature order.

## Software and one command

Use Python 3.11, R (reference versions in `validation/versiones_paquetes.csv` and `validation/sesion_R.txt`), CmdStan 2.36.0, a C++17 compiler and GNU make. PDF generation additionally needs TeX Live with `latexmk`. Install the R dependencies once:

```r
install.packages(c("sf", "sfnetworks", "igraph", "Matrix", "RSpectra", "jsonlite",
                  "posterior", "loo", "statmod", "ggplot2", "dplyr", "tidyr", "scales"))
install.packages("cmdstanr", repos=c("https://stan-dev.r-universe.dev", getOption("repos")))
cmdstanr::install_cmdstan(version="2.36.0")
```

From this directory, run:

```sh
python3 reproduce.py --output run --cmdstan "$HOME/.cmdstan/cmdstan-2.36.0" --workers 4
```

`--output` must name a new directory. The installation location is the only machine-specific input; project paths are derived from the delivered files. Run from any working directory using the path to `reproduce.py`. `--skip-pdf` omits TeX compilation while still generating every numerical result, table and figure. Four workers run the four chains of one fit concurrently; 20 fits run in sequence. Allow several GB of disk and RAM for chain output and posterior summaries. Runtime depends on the machine; the four BYM2 fits are substantially slower.

The reference settings fix seeds, initializations, four chains, warmup and sampling budgets (1,000 + 1,000 per chain for other fits; 2,000 + 2,000 for BYM2), `adapt_delta` (0.99 or 0.995), maximum depth (13 or 14) and model parameters. There are no GUI analysis steps and no undocumented input parameters. Stan compilation is local; no platform-specific binary is required. Reference results are never treated as newly fitted results.

## Complete data workflow

1. Verify checksums, load record projections and fixed road snapshots. Validate coordinates and check that repeated incident identifiers do not disagree on location. Do **not** deduplicate the record-count response.
2. Barcelona: transform WGS84 points to EPSG:25831, retain the analysis-road bounding box, assign each record to the nearest of 2,714 fixed road features and aggregate.
3. Bogotá: transform to EPSG:3116, retain the full source-road bounding box, assign each row to its nearest feature among 54,497 frozen source roads, aggregate by OSM ID and transfer counts to the 8,169 selected road features. This reproduces the original pre-clipping aggregation, including contributions on clipped-away portions of an inherited OSM feature; it is not a new within-polygon crash-incidence definition. Assignment IDs and distances are exported for audit.
4. Create the undirected network from fixed geometries; retain the component with the greatest number of edges, preserve parallel features and source order. Construct binary edge adjacency from shared endpoints (a pair sharing both endpoints receives one adjacency). Retained supports have 2,540 and 6,079 segments, with 3,735 and 85,068 records.
5. Compute lengths in km and eight attribute contrasts. Treat invalid or missing lane/surface tags as documented unknown categories. Compute 61 leading normalized-adjacency eigenpairs, fix signs by largest absolute entry, omit the Perron vector and center/RMS-scale the next 20 candidates. Verify eigen-equation residuals. Save graph, response, design and basis.
6. Generate 20 model JSON inputs; apply the deposited final settings and initializations; compile the shared Stan model and sample every fit. No results are reconstructed from typed table numbers.
7. Summarize actual draws: coefficients and ratios, slab membership and latent K distribution, signed contributions and fields; sampling diagnostics; posterior predictive checks; residual Moran statistics; pointwise PSIS-LOO. BYM2 LOO integrates independent effects using checked Gaussian quadrature. Paired ELPD differences and their standard errors use segmentwise differences.
8. Export all metrics (`run/metrics/reported_findings.csv`), all nine table bodies, six figure PDFs and the compiled document PDFs. `run/verification.json` compares new numerical summaries with references; floating point/compiler/library differences can prevent bitwise equality. Sampling warnings are retained, not automatically declared acceptable.

Outputs include `run/source_reconstruction/audits` (record assignments and exclusions), `run/modelos_espectrales/datos` (derived model inputs), `run/work/modelos_espectrales/cadenas` (actual chains), `run/modelos_espectrales/resultados` (fit summaries), and `run/documents` (tables, figures and PDFs).

## Finding-by-finding instructions

Run the command above once; it executes source processing and inference for **every row below**. The final column identifies the draw-derived result used by the exporter. `C` below means Barcelona or Bogota; specification suffixes are `sin_covariables` (U), `con_covariables` (C), and `con_covariables__longitud_libre` (L). Models are `NB`, `CAR_espectral` (dense), `Seleccion_espectral` (mixture) and `BYM2_espectral`. L is fitted for NB and mixture only. Fit directory names join city, model and specification with `__`.

| Finding | Generated artifact under `run/documents` | Derived numerical inputs under `run/modelos_espectrales` |
|---|---|---|
| Table 1: segments, intersections, total records, lengths and zeros | `tables/data.tex` | `resultados/auditoria_datos.csv`, derived from reconstructed records and retained graph |
| Table 2: expected K, 90% interval and membership above 0.5 | `tables/membership.tex` | mixture fits: `seleccion_modos.csv`, `distribucion_K.csv`; K quantiles use the latent membership distribution |
| Table 3: count ratios and intervals | `tables/effects.tex` | mixture L: `covariables.csv`, posterior coefficient transformations |
| Table 4: ELPD, uncertainty, Pareto and sampler flags | `tables/prediction.tex` | all 20 fits: `prediccion_loo.csv`, `diagnosticos.csv`, `loo.rds` |
| Figure 1: eigenvalue spectra | `figures/spectrum.pdf` | `resultados/C_espectro.csv` from retained adjacency |
| Figure 2: first three geographical candidate modes | `figures/raw_network_modes.pdf` | `datos/C_preparado.rds`: raw eigenvectors and retained geometry |
| Figure 3: mode membership U/C/L | `figures/mode_membership.pdf` | six mixture fits: `seleccion_modos.csv`; aggregate `seleccion_con_sensibilidad_longitud.csv` |
| Figure 4: three greatest RMS signed contributions | `figures/adjusted_contributions.pdf` | mixture L: `contribuciones_modos.csv` and retained geometry/basis |
| Figure 5: full posterior mean field | `figures/adjusted_field.pdf` | mixture L: `resumen_segmentos.csv` and retained geometry |
| Table S1: recorded attribute categories | `tables/s1.tex` | `resultados/covariables_niveles.csv` from fixed road tags and retained graph |
| Table S2: all mixture L count ratios | `tables/s2.tex` | mixture L: `covariables.csv`, including nuisance unknown categories |
| Table S3: core diagnostics | `tables/s3.tex` | all fits: `diagnosticos.csv`, actual chains |
| Table S4: Bayesian predictive checks | `tables/s4.tex` | mixture L: `comprobaciones_predictivas.csv`, actual posterior predictive draws |
| Table S5: residual Moran checks | `tables/s5.tex` | mixture L: `autocorrelacion_residual.csv`, adjacency and posterior predictive draws |
| Figure S1: marginal field interval signs | `figures/field_uncertainty.pdf` | mixture L: `resumen_segmentos.csv`, 5th/95th field percentiles |
| Narrative length elasticities and intervals | `run/metrics/reported_findings.csv` | four L fits: `elasticidad_longitud.csv` |
| Narrative covariate and length predictive gains; mixture versus NB L | same metrics file, `Narrative paired ELPD` | segmentwise differences of the corresponding `loo.rds` pointwise vectors; SE = sqrt(n × variance(difference)) |
| Narrative selected modes, K, dominant contributions, ratios, residual checks and posterior predictive values | same metrics file, labelled by fit, metric and output type | all mode membership/MCSE/ESS/Rhat rows; K distribution; coefficient summaries; contribution RMS; predictive moments/zeros/tails; Moran statistics |
| Numerical eigenpair checks and field interval counts | same metrics file | all 61 eigenvalues/residuals; mixture L posterior field percentiles |

Every table and figure in the current supplementary PDF is included. It is limited to three pages. The older 1,200-dataset/4,800-fit simulation was not recovered and is **not a reported experiment or supporting result of the current manuscript**; archived earlier material must not be mistaken for current findings. No historical numerical simulation result is manufactured by this workflow.

## Execution evidence and limits

`validation/` records a genuine rerun of all 20 models, source-to-model input comparisons and a collector of 4,104 numerical assessment rows. Comparison of 172 CSV summaries checked 613,334 numerical cells with zero differences outside relative tolerance 1e-7 / absolute tolerance 1e-8. Separate reconstruction from the original rows matched the response, offset, design matrix, candidate bases, eigenvalues, adjacency and endpoints exactly in both cities. New model results were sampled in an empty output workspace. The final export used the new summaries. This verifies the computational reconstruction; it does not resolve the four BYM2 energy warnings, influential PSIS cases, residual dependence, missing traffic volume, unequal source units/periods or untested geographical transfer. These limitations remain in the manuscript.
