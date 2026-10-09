import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from .config import ROOT, STAGES

def save(name):
    path = ROOT / 'figures' / (name + '.png')
    path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(path, dpi=160, bbox_inches='tight')
    plt.close('all')

def embeddings(a, name, keys, basis='umap'):
    for key in keys:
        sc.pl.embedding(a, basis=basis, color=key, show=False, frameon=False,
                        title=f'{name.replace("_", " ")} | {key}', size=3)
        save(f'umap/{name}_{key}')

def qc_plot(frame):
    fig, axs = plt.subplots(1, 3, figsize=(16, 4))
    samples = list(frame.sample_id.unique())
    for ax, key in zip(axs, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt']):
        ax.boxplot([frame.loc[frame.sample_id == s, key] for s in samples], showfliers=False)
        ax.set_xticks(range(1, len(samples) + 1), [s[-3:] for s in samples], rotation=45)
        ax.set_xlabel('GSM accession suffix'); ax.set_ylabel(key)
        if key == 'n_genes_by_counts':
            ax.axhline(200, color='red', ls='--'); ax.axhline(5000, color='red', ls='--')
        elif key == 'pct_counts_mt':
            ax.axhline(15, color='red', ls='--')
        else:
            ax.set_yscale('log')
    fig.suptitle('Before project QC: supplied GEO matrix cells, per sample')
    save('qc/per_sample')
    fig, axs = plt.subplots(3, 3, figsize=(12, 10))
    for ax, s in zip(axs.flat, samples):
        d = frame[frame.sample_id == s]
        ax.scatter(d.total_counts, d.n_genes_by_counts, c=d.pct_counts_mt,
                   s=1, vmin=0, vmax=30, cmap='viridis', rasterized=True)
        ax.set_xscale('log'); ax.set_title(s); ax.set_xlabel('Counts'); ax.set_ylabel('Genes')
    fig.tight_layout(); save('qc/counts_vs_genes')

def composition(obs, category, name):
    counts = pd.crosstab(obs.sample_id, obs[category])
    fractions = counts.div(counts.sum(axis=1), axis=0)
    meta = obs[['sample_id', 'disease_stage']].drop_duplicates().set_index('sample_id')
    fractions.to_csv(ROOT / 'results' / f'{name}_sample_fractions.csv')
    counts.to_csv(ROOT / 'results' / f'{name}_sample_counts.csv')
    labels = [f'{s[-3:]} | {meta.loc[s, "disease_stage"]}' for s in fractions.index]
    fractions.plot.bar(stacked=True, figsize=(12, 5), colormap='tab20')
    plt.xticks(range(len(labels)), labels, rotation=55, ha='right')
    plt.ylabel('Fraction within sample'); plt.xlabel('Sample'); plt.ylim(0, 1)
    plt.legend(bbox_to_anchor=(1.02, 1), title=category)
    save(f'progression/{name}_by_sample')
    # Each sample contributes equally; stage totals must not weight big samples more.
    stage = fractions.join(meta).groupby('disease_stage', observed=True).mean().reindex(STAGES)
    stage.to_csv(ROOT / 'results' / f'{name}_stage_mean_fractions.csv')
    stage.plot.bar(stacked=True, figsize=(10, 5), colormap='tab20')
    plt.ylabel('Mean of within-sample fractions'); plt.xlabel('Disease group')
    plt.xticks(rotation=20); plt.ylim(0, 1); plt.legend(bbox_to_anchor=(1.02, 1), title=category)
    save(f'progression/{name}_by_stage')
