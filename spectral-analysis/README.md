# Bayesian spectral analysis on road networks

Use the [finding-by-finding reproduction guide](REPRODUCIBILITY.md), then run the delivered [source-to-results executor](reproduction/reproduce.py). The self-contained package is in `reproduction/`; it includes source-record projections, fixed road snapshots, source hashes, preprocessing, Stan code, recorded settings, initializations, numerical references and document sources. [Download everything](https://github.com/DannaCruz/CrashBayesian/archive/refs/heads/main.zip).

The current documents are the [main analysis](documents/manuscript.pdf) and [three-page supplement](documents/supplement.pdf). All four main tables, five main figures, five supplementary tables, the supplementary figure and narrative assessment metrics are mapped to their computational outputs. Historical simulations are not reported findings of these documents.

The 20 newly executed fits matched 613,334 numerical cells from 172 reference summaries within the stated tolerance. See [execution evidence](reproduction/validation/numeric_comparison.json) and [source-to-model identity checks](reproduction/validation/source_to_model_input_identity.csv). Sixteen fits passed the joint sampler criterion; four BYM2 fits retain diagnostic warnings. Computational reproduction does not establish model adequacy or geographical transfer.

Barcelona uses accident-report rows; Bogotá uses person-involvement rows. These are different responses and periods, not comparable counts of distinct crashes or traffic-adjusted risks. The guide describes the preserved pre-clipping aggregation in Bogotá and remaining substantive limitations.

`analysis/` retains the earlier aggregate-input runner and reference outputs for inspection. The complete supported workflow is now `reproduction/reproduce.py`; it replaces external-network-file requirements in the older scripts. The separate node–edge package uses different data and specifications.
