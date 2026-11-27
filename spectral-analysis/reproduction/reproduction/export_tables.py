"""Export empirical LaTeX table bodies from saved numerical summaries.

This command does not fit models, reconstruct likelihoods, or validate original
municipal records. Use a fresh output directory; reference data are never changed.
"""
from pathlib import Path
import argparse
import csv

HERE = Path(__file__).resolve().parent
CITIES = ['Barcelona', 'Bogota']
MODELS = ['NB', 'CAR_espectral', 'Seleccion_espectral', 'BYM2_espectral']
SPECS = ['sin_covariables', 'con_covariables', 'con_covariables__longitud_libre']
LABELS = dict(zip(SPECS, ['U', 'C', 'L']))
COVARIATES = {
    'highwaytertiary': 'Tertiary road',
    'highwaysecondary': 'Secondary road',
    'highwayprimary': 'Primary road',
    'lanesone': 'One lane',
    'lanesthree_plus': 'Three or more lanes',
    'lanesunknown': 'Unknown lane count',
    'surfaceother_recorded': 'Other surface tag',
    'surfaceunknown': 'Unknown surface',
}

def read(path):
    with path.open(newline='', encoding='utf-8') as stream:
        return list(csv.DictReader(stream))

def city_label(city):
    return r"Bogot\'a" if city == 'Bogota' else city

def fmt(value, digits=2):
    return f'{float(value):.{digits}f}'

def truth(value):
    return str(value).lower() in ['true', '1']

def build(fits):
    def fit(city, model='Seleccion_espectral', spec=SPECS[2]):
        return fits / f'{city}__{model}__{spec}'
    def one(city, model, spec, filename):
        rows = read(fit(city, model, spec) / filename)
        if len(rows) != 1:
            raise ValueError(f'Expected one row: {city}/{model}/{spec}/{filename}')
        return rows[0]

    tables = {}
    audits = {row['city']: row for row in read(fits / 'auditoria_datos.csv')}
    rows = []
    for label, key, digits in [
        ('Retained segments', 'retained_segments', None),
        ('Retained intersections', 'retained_intersections', None),
        ('Recovered total count', 'retained_count_total', None),
        ('Represented length (km)', 'length_km', 1),
    ]:
        values = [f"{int(audits[c][key]):,}" if digits is None else fmt(audits[c][key], digits) for c in CITIES]
        rows.append([label, *values])
    rows.append([r'Zero-count segments (\%)', *[fmt(100 * int(audits[c]['zero_segments']) / int(audits[c]['retained_segments']), 1) for c in CITIES]])
    tables['data.tex'] = rows

    rows = []
    for city in CITIES:
        for spec in SPECS:
            modes = read(fit(city, spec=spec) / 'seleccion_modos.csv')
            pk = sorted(read(fit(city, spec=spec) / 'distribucion_K.csv'), key=lambda x: int(x['K']))
            assert abs(sum(float(x['probability']) for x in pk) - 1) < 1e-7
            # Quantiles of latent K use its posterior mixture distribution, not the
            # draw-wise sum of membership probabilities.
            def quantile(p):
                total = 0.0
                for row in pk:
                    total += float(row['probability'])
                    if total >= p:
                        return int(row['K'])
                raise ValueError('Invalid K distribution')
            active = ', '.join(str(int(x['mode'])) for x in modes if float(x['PIP']) > .5) or 'None'
            rows.append([city_label(city), LABELS[spec], fmt(sum(float(x['PIP']) for x in modes)), f'[{quantile(.05)}, {quantile(.95)}]', active])
    tables['membership.tex'] = rows

    effects = {c: {x['covariate']: x for x in read(fit(c) / 'covariables.csv')} for c in CITIES}
    substantive = ['highwaytertiary','highwaysecondary','highwayprimary','lanesone','lanesthree_plus','surfaceother_recorded']
    contrasts = ['Tertiary vs residential road','Secondary vs residential road','Primary vs residential road','One vs two lanes','Three or more vs two lanes','Other surface tag vs asphalt']
    rows = []
    for key, label in zip(substantive, contrasts):
        values = []
        for c in CITIES:
            x = effects[c][key]
            values.append(f"{fmt(x['rate_ratio_median'])} [{fmt(x['rate_ratio_q05'])}, {fmt(x['rate_ratio_q95'])}]")
        rows.append([label, *values])
    tables['effects.tex'] = rows

    prediction = []
    diagnostics = []
    model_labels = ['NB', 'Dense spectral', 'Mixture spectral', 'Spectral BYM2']
    for c in CITIES:
        prediction.append([r'\textit{' + city_label(c) + '}', '', '', '', '', ''])
        city_diagnostics = []
        for typ, model in enumerate(MODELS):
            for spec in SPECS:
                if spec == SPECS[2] and model not in ['NB','Seleccion_espectral']:
                    continue
                pr = one(c, model, spec, 'prediccion_loo.csv')
                dg = one(c, model, spec, 'diagnosticos.csv')
                status = 'Pass' if truth(dg['diagnostic_pass']) else 'Flagged'
                prediction.append([model_labels[typ], LABELS[spec], fmt(pr['elpd_loo'],1), fmt(pr['se_elpd'],1), int(pr['pareto_bad']), status])
                city_diagnostics.append((f'{c}__{model}__{spec}', [city_label(c), ['N','D','M','B'][typ] + '-' + LABELS[spec], dg['draws_per_chain'], fmt(dg['max_rhat'],3), fmt(dg['min_bulk_ess'],0), fmt(dg['min_ebfmi'],3), status]))
        # Match the archived diagnostic table's lexicographic fit-ID ordering.
        diagnostics.extend(row for _, row in sorted(city_diagnostics))
        prediction.append([r'\addlinespace'])
    tables['prediction.tex'] = prediction
    assert len(diagnostics) == 20
    tables['s3.tex'] = diagnostics

    categories = read(fits / 'covariables_niveles.csv')
    rows = []
    for var, levels, label in [
        ('highway', ['residential','tertiary','secondary','primary'], 'Road class'),
        ('lanes', ['two','one','three_plus','unknown'], 'Lanes'),
        ('surface', ['asphalt','other_recorded','unknown'], 'Surface tag'),
    ]:
        for level in levels:
            values = [next(x['n'] for x in categories if x['city']==c and x['variable']==var and x['level']==level) for c in CITIES]
            rows.append([label, level.replace('_',r'\_'), *values])
    tables['s1.tex'] = rows
    tables['s2.tex'] = [[city_label(c), label, *[fmt(effects[c][key][k]) for k in ['rate_ratio_median','rate_ratio_q05','rate_ratio_q95']]] for c in CITIES for key,label in COVARIATES.items()]

    rows = []
    stat_labels = dict(zero_fraction='Zero fraction',variance='Variance',maximum='Maximum',q99='99th percentile',mean='Mean')
    for c in CITIES:
        for x in read(fit(c) / 'comprobaciones_predictivas.csv'):
            n = 3 if x['statistic'] in ['zero_fraction','mean'] else 2
            rows.append([city_label(c), stat_labels[x['statistic']], fmt(x['observed'],n), fmt(x['replicated_median'],n), f"[{fmt(x['replicated_q05'],n)}, {fmt(x['replicated_q95'],n)}]", fmt(x['bayesian_p_upper'],3)])
    tables['s4.tex'] = rows
    rows = []
    for c in CITIES:
        for model, label in [('NB','N-L'),('Seleccion_espectral','M-L')]:
            x = one(c, model, SPECS[2], 'autocorrelacion_residual.csv')
            rows.append([city_label(c), label, fmt(x['observed_I_median'],3), fmt(x['replicated_I_median'],3), f"[{fmt(x['replicated_I_q05'],3)}, {fmt(x['replicated_I_q95'],3)}]", fmt(x['bayesian_p_upper'])])
    tables['s5.tex'] = rows
    return tables

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fits', type=Path, default=HERE / 'reference/fits')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Output directory exists; choose a new directory.')
    # Read and validate all required summaries before creating any output.
    tables = build(args.fits)
    args.output.mkdir(parents=True)
    for name, rows in tables.items():
        lines = []
        for row in rows:
            if row == [r'\addlinespace']:
                lines.append(row[0])
            else:
                lines.append(' & '.join(map(str,row)) + r' \\')
        (args.output / name).write_text('\n'.join(lines) + '\n\\bottomrule\n', encoding='utf-8')
    (args.output / 'SCOPE.txt').write_text('Exported nine empirical table bodies from saved numerical CSV summaries. No models were fitted, likelihoods recomputed, original events verified or historical simulations reproduced.\n')
    print('Exported Tables 1–4 and S1–S5. No sampling performed.')

if __name__ == '__main__':
    main()
