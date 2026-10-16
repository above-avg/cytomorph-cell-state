"""Sample-level host-response contrasts and an audit of deposited viral features.

Run with python -m src.host_response after the primary pipeline.
No infection status is inferred for individual cells and no integration sites are called.
"""
import gzip
import hashlib
import json
import re
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

from .config import ROOT, SIGNATURES, STAGES

PROGRAMS = ['Interferon_response', 'Cell_cycle', 'Cytotoxicity', 'Inhibitory_receptors']
CONTRASTS = {
    'HPV_status_in_normal': ([STAGES[1]], [STAGES[0]]),
    'HSIL_vs_HPV_positive_normal': ([STAGES[2]], [STAGES[1]]),
    'Cancer_vs_HSIL': ([STAGES[3]], [STAGES[2]]),
    'Cancer_vs_noncancer': ([STAGES[3]], STAGES[:3]),
}
MIN_CELLS = 20


def file_digest(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def contrast(values, target, reference):
    """Equal sample weights; omission is a sensitivity range, not a confidence interval."""
    if set(target) & set(reference):
        raise ValueError('Contrast groups overlap')
    if len(set(target + reference)) != len(target + reference):
        raise ValueError('Repeated sample IDs')
    if len(target) < 2 or len(reference) < 2:
        raise ValueError('At least two samples in each original group are required')
    if not values.index.is_unique or not np.isfinite(values.loc[target + reference]).all():
        raise ValueError('Expected unique samples with finite scores')
    full = float(values.loc[target].mean() - values.loc[reference].mean())
    folds = []
    for omitted in target + reference:
        t = [s for s in target if s != omitted]
        r = [s for s in reference if s != omitted]
        effect = float(values.loc[t].mean() - values.loc[r].mean())
        folds.append(dict(omitted_sample=omitted, effect=effect,
                          n_target=len(t), n_reference=len(r)))
    effects = np.array([f['effect'] for f in folds])
    return dict(effect=full, omission_min=float(effects.min()), omission_max=float(effects.max()),
                direction_retained=bool(full != 0 and (np.sign(effects) == np.sign(full)).all()),
                max_absolute_change=float(np.abs(effects - full).max())), folds


def viral_feature_matches(features):
    # Explicit virus labels only: E6/E7 alone can collide with host feature names.
    pattern = r'(?i)(?:(?:^|[^a-z0-9])HPV(?:[-_ ]?\d+)?(?=[^a-z0-9]|$)|human[_ -]*papillomavirus|papillomavirus)'
    return features.astype(str).apply(lambda col: col.str.contains(pattern, regex=True)).any(axis=1)


def main():
    out = ROOT / 'results/host_response'
    figdir = ROOT / 'figures/host_response'
    out.mkdir(parents=True, exist_ok=True)
    figdir.mkdir(parents=True, exist_ok=True)
    samples = pd.read_csv(ROOT / 'data/metadata/samples.csv').set_index('sample_id')
    checkpoint = ROOT / 'data/processed/annotated.h5ad'
    print('Loading annotated checkpoint', flush=True)
    a = sc.read_h5ad(checkpoint)
    coverage = []
    obs = a.obs[['sample_id', 'cell_type']].copy()
    for program in PROGRAMS:
        present = [g for g in SIGNATURES[program] if g in a.var_names]
        missing = [g for g in SIGNATURES[program] if g not in a.var_names]
        if not present:
            raise ValueError(f'No genes available for {program}')
        obs[program] = np.asarray(a[:, present].X.mean(axis=1)).ravel()
        coverage.append(dict(program=program, present=';'.join(present), missing=';'.join(missing)))
    summary = obs.groupby(['sample_id', 'cell_type'], observed=True)[PROGRAMS].mean()
    summary['n_cells'] = obs.groupby(['sample_id', 'cell_type'], observed=True).size()
    summary = summary.reset_index().merge(samples.reset_index(), on='sample_id', validate='many_to_one')
    summary['eligible'] = (summary.n_cells >= MIN_CELLS) & (summary.cell_type != 'Ambiguous')
    summary.to_csv(out / 'sample_cell_type_scores.csv', index=False)
    pd.DataFrame(coverage).to_csv(out / 'gene_coverage.csv', index=False)
    print('Computing sample omission contrasts', flush=True)
    results, folds, skipped = [], [], []
    for cell_type, group in summary[summary.eligible].groupby('cell_type', observed=True):
        group = group.set_index('sample_id')
        for name, (target_stages, reference_stages) in CONTRASTS.items():
            target = list(group.index[group.disease_stage.isin(target_stages)])
            reference = list(group.index[group.disease_stage.isin(reference_stages)])
            for program in PROGRAMS:
                keys = dict(cell_type=str(cell_type), contrast=name, program=program)
                if len(target) < 2 or len(reference) < 2:
                    skipped.append(dict(**keys, n_target=len(target), n_reference=len(reference),
                                        reason='Fewer than two eligible samples in a group'))
                    continue
                result, omitted = contrast(group[program], target, reference)
                # Excluding ADC is a histology sensitivity, not a separate validation cohort.
                scc = [s for s in target if group.loc[s, 'filename_histology'] == 'SCC']
                scc_effect = float(group.loc[scc, program].mean() - group.loc[reference, program].mean()) if target_stages == ['Cancer'] and len(scc) >= 2 else np.nan
                results.append(dict(**keys, n_target=len(target), n_reference=len(reference),
                                    **result, SCC_only_effect=scc_effect))
                folds.extend(dict(**keys, **f) for f in omitted)
    effects = pd.DataFrame(results)
    effects.to_csv(out / 'contrasts.csv', index=False)
    pd.DataFrame(folds).to_csv(out / 'sample_omission_folds.csv', index=False)
    pd.DataFrame(skipped, columns=['cell_type', 'contrast', 'program', 'n_target', 'n_reference', 'reason']).to_csv(out / 'skipped_contrasts.csv', index=False)

    # Audit original, pre-QC feature files rather than only retained host genes.
    audit = []
    for sample, row in samples.iterrows():
        path = ROOT / 'data/raw' / (row['prefix'] + '.features.tsv.gz')
        with gzip.open(path, 'rt') as stream:
            features = pd.read_csv(stream, sep='\t', header=None)
        matches = features[viral_feature_matches(features)]
        audit.append(dict(sample_id=sample, features=len(features), explicitly_labelled_viral_features=len(matches),
                          matching_features=';'.join('|'.join(map(str, r)) for r in matches.itertuples(index=False, name=None)),
                          feature_file_sha256=file_digest(path)))
    pd.DataFrame(audit).to_csv(out / 'viral_feature_audit.csv', index=False)

    # Direct normal HPV-status comparison and subsequent lesion contrasts side by side.
    fig, axes = plt.subplots(1, 4, figsize=(17, 5), sharey=True, sharex=True)
    types = sorted(effects.cell_type.unique())
    for ax, name in zip(axes, CONTRASTS):
        d = effects[(effects.contrast == name) & (effects.program == 'Interferon_response')].set_index('cell_type')
        for y, label in enumerate(types):
            if label not in d.index:
                continue
            row = d.loc[label]
            ax.plot([row.omission_min, row.omission_max], [y, y], color='#6689a8', lw=3)
            ax.scatter(row.effect, y, color='#102f4d', zorder=3)
            if pd.notna(row.SCC_only_effect):
                ax.scatter(row.SCC_only_effect, y, marker='x', color='#bd590c', zorder=4)
        ax.axvline(0, color='gray', ls='--', lw=1)
        ax.set_title(name.replace('_', ' '), fontsize=10)
        ax.set_xlabel('Difference in sample-mean score')
        ax.set_yticks(range(len(types)), types)
    fig.suptitle('Host interferon response: dots = full contrast; lines = sample-omission range; orange x = SCC only')
    fig.tight_layout()
    fig.savefig(figdir / 'interferon_sample_influence.png', dpi=160, bbox_inches='tight')
    plt.close(fig)

    selected = summary[(summary.cell_type == 'T/NK') & summary.eligible]
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    for ax, program in zip(axes, ['Interferon_response', 'Inhibitory_receptors']):
        for x, stage in enumerate(STAGES):
            d = selected[selected.disease_stage == stage].sort_values('sample_id')
            for jitter, (_, row) in zip(np.linspace(-.16, .16, len(d)), d.iterrows()):
                ax.scatter(x + jitter, row[program], color='#236b83', s=45)
                ax.annotate(row.sample_id[-3:], (x + jitter, row[program]), xytext=(4, 4), textcoords='offset points', fontsize=8)
        ax.set_xticks(range(4), ['Normal HPV-', 'Normal HPV+', 'HSIL', 'Cancer'], rotation=15)
        ax.set_title(program.replace('_', ' '))
        ax.set_ylabel('Mean log-normalized expression')
        ax.set_xlim(-.5, 3.55)
        ax.margins(y=.2)
    fig.suptitle('T/NK host programs: one point per sample; labels are GSM accession suffixes')
    fig.tight_layout()
    fig.savefig(figdir / 'T_NK_sample_programs.png', dpi=160, bbox_inches='tight')
    plt.close(fig)
    manifest = dict(status='completed', samples=len(samples), min_cells_per_sample_cell_type=MIN_CELLS,
                    evaluated_contrasts=len(effects), skipped_contrasts=len(skipped),
                    omission_folds=len(folds), eligible_sample_cell_type_rows=int(summary.eligible.sum()),
                    cells=a.n_obs, programs=PROGRAMS, checkpoint_sha256=file_digest(checkpoint),
                    sample_metadata_sha256=file_digest(ROOT / 'data/metadata/samples.csv'),
                    interpretation='Conditional sample omission; no refitting of annotation or normalization; no causal HPV or integration claim')
    (out / 'manifest.json').write_text(json.dumps(manifest, indent=2))
    print(json.dumps(manifest, indent=2), flush=True)


if __name__ == '__main__':
    main()
