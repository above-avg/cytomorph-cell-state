import json
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from sklearn.metrics import adjusted_rand_score
from .config import ROOT, STAGES, MARKERS, SIGNATURES, SEED
from .plotting import save, embeddings, composition
from .preprocessing import graph

def annotate(a):
    available = {k: [g for g in v if g in a.var_names] for k, v in MARKERS.items()}
    for label, genes in available.items():
        if len(genes) < 2:
            raise ValueError(f'Insufficient marker coverage for {label}')
        sc.tl.score_genes(a, genes, score_name='broad_' + label, random_state=SEED, use_raw=True)
    keys = ['broad_' + k for k in available]
    scores = a.obs.groupby('leiden', observed=True)[keys].mean()
    mapping, review = {}, []
    for cluster, values in scores.iterrows():
        ranked = values.sort_values(ascending=False)
        label = ranked.index[0].removeprefix('broad_')
        margin = float(ranked.iloc[0] - ranked.iloc[1])
        if ranked.iloc[0] <= 0 or margin < 0.1:
            label = 'Ambiguous'
        mapping[cluster] = label
        review.append(dict(cluster=cluster, provisional_label=label, top_score=float(ranked.iloc[0]),
                           margin=margin, second_label=ranked.index[1].removeprefix('broad_'), cells=int((a.obs.leiden == cluster).sum())))
    a.obs['cell_type'] = pd.Categorical(a.obs.leiden.astype(str).map(mapping))
    scores.to_csv(ROOT / 'results/annotation_scores.csv')
    pd.DataFrame(review).to_csv(ROOT / 'results/annotation_review.csv', index=False)
    sc.pl.dotplot(a, available, groupby='leiden', use_raw=True, standard_scale='var', show=False)
    save('annotation/broad_marker_dotplot')
    composition(a.obs, 'cell_type', 'broad_composition')
    embeddings(a, 'all_cells', ['cell_type', 'disease_stage', 'sample_id', 'leiden'])
    embeddings(a, 'all_cells_baseline', ['cell_type', 'disease_stage', 'sample_id'], basis='umap_uncorrected')
    return a

def markers(a, group, prefix):
    print('Marker ranking', prefix, flush=True)
    # Descriptive cell-level rankings, NOT patient-level disease significance.
    # Restrict ranking to the 500 most variable genes. Ranking every detected
    # gene across ~69k cells is needlessly expensive; this prespecified subset
    # keeps candidate marker discovery tractable and reproducible.
    work = a[:, a.var['highly_variable']].copy()
    keep = np.argsort(work.var['dispersions_norm'].to_numpy())[-500:]
    work = work[:, keep].copy()
    sc.tl.rank_genes_groups(work, groupby=group, method='t-test', use_raw=False, pts=True)
    frame = sc.get.rank_genes_groups_df(work, group=None)
    frame.to_csv(ROOT / f'results/markers/{prefix}_all.csv.gz', index=False)
    frame.groupby('group', observed=True).head(20).to_csv(ROOT / f'results/markers/{prefix}_top20.csv', index=False)

def signatures(a, prefix):
    audit = []
    keys = []
    for label, genes in SIGNATURES.items():
        present = [g for g in genes if g in a.var_names]
        audit.append(dict(signature=label, requested=';'.join(genes), present=';'.join(present), missing=';'.join(g for g in genes if g not in present)))
        if present:
            key = 'signature_' + label
            # Transparent arithmetic mean of log-normalized expression, not enrichment.
            a.obs[key] = np.asarray(a[:, present].X.mean(axis=1)).ravel()
            keys.append(key)
    pd.DataFrame(audit).to_csv(ROOT / f'results/signatures/{prefix}_gene_coverage.csv', index=False)
    sample = a.obs.groupby(['sample_id', 'disease_stage'], observed=True)[keys].mean().reset_index()
    sample.to_csv(ROOT / f'results/signatures/{prefix}_sample_means.csv', index=False)
    fig, axs = plt.subplots(2, 4, figsize=(17, 8))
    for ax, key in zip(axs.flat, keys):
        for stage_id, stage in enumerate(STAGES):
            d = sample[sample.disease_stage == stage]
            ax.scatter(np.full(len(d), stage_id) + np.linspace(-.12, .12, len(d)), d[key], s=45)
        ax.set_xticks(range(4), ['N-', 'N+', 'HSIL', 'CA']); ax.set_title(key.removeprefix('signature_'))
        ax.set_ylabel('Sample mean log-normalized expression')
    fig.suptitle(f'{prefix}: each point is one sample; signatures prespecified, exploratory')
    fig.tight_layout(); save(f'progression/{prefix}_signature_scores')
    return sample

def integration_diagnostics(a, prefix):
    # Compare membership, rather than equating a visually mixed UMAP with success.
    report = dict(cells=a.n_obs, baseline_clusters=int(a.obs.leiden_uncorrected.nunique()),
                  harmony_clusters=int(a.obs.leiden.nunique()),
                  adjusted_rand_index=float(adjusted_rand_score(a.obs.leiden_uncorrected, a.obs.leiden)))
    for key in ['leiden', 'leiden_uncorrected']:
        counts = pd.crosstab(a.obs[key], a.obs.sample_id)
        counts.to_csv(ROOT / f'results/{prefix}_{key}_sample_counts.csv')
        report[key + '_median_dominant_sample_fraction'] = float(counts.div(counts.sum(axis=1), axis=0).max(axis=1).median())
    (ROOT / f'results/{prefix}_integration_diagnostics.json').write_text(json.dumps(report, indent=2))

def epithelial(a):
    e = a[a.obs.cell_type == 'Epithelial'].copy()
    if e.n_obs < 100:
        raise ValueError('Too few provisional epithelial cells; revise annotation')
    # Recompute HVGs/graph within the subset, retaining full expression for markers.
    sc.pp.highly_variable_genes(e, n_top_genes=2000, flavor='seurat', batch_key='sample_id')
    graph(e, 'epithelial', resolution=0.5)
    embeddings(e, 'epithelial', ['leiden', 'disease_stage', 'sample_id'])
    embeddings(e, 'epithelial_baseline', ['disease_stage', 'sample_id'], basis='umap_uncorrected')
    composition(e.obs, 'leiden', 'epithelial_composition')
    integration_diagnostics(e, 'epithelial')
    markers(e, 'leiden', 'epithelial')
    signatures(e, 'epithelial')
    e.write_h5ad(ROOT / 'data/processed/epithelial.h5ad', compression='gzip')
    return e

def immune(a):
    subset = a[a.obs.cell_type.isin(['T/NK', 'B', 'Plasma', 'Myeloid', 'Mast'])].copy()
    signatures(subset, 'immune')
    # T/NK-specific receptor/cytotoxicity avoids confounding with other immune types.
    t = subset[subset.obs.cell_type == 'T/NK'].copy()
    if t.n_obs:
        signatures(t, 'T_NK')
    composition(subset.obs, 'cell_type', 'immune_composition')

def pseudobulk_and_lopo(e, samples):
    print('Pseudobulk and held-out sample analysis', flush=True)
    vectors, info = [], []
    genes = list(e.var_names)
    for row in samples.to_dict('records'):
        ids = e.obs_names[e.obs.sample_id == row['sample_id']]
        if len(ids) < 20:
            continue
        raw = sc.read_h5ad(ROOT / f'data/processed/{row["sample_id"]}_counts.h5ad')
        v = np.asarray(raw[ids, genes].X.sum(axis=0)).ravel().astype('float64')
        vectors.append(v)
        info.append(dict(sample_id=row['sample_id'], disease_stage=row['disease_stage'], cells=len(ids)))
    counts = pd.DataFrame(vectors, index=[r['sample_id'] for r in info], columns=genes)
    meta = pd.DataFrame(info).set_index('sample_id')
    counts.to_csv(ROOT / 'results/signatures/epithelial_pseudobulk_counts.csv.gz')
    meta.to_csv(ROOT / 'results/signatures/epithelial_pseudobulk_samples.csv')
    logcpm = np.log1p(counts.div(counts.sum(axis=1), axis=0) * 1e6)
    logcpm.to_csv(ROOT / 'results/signatures/epithelial_pseudobulk_logcpm.csv.gz')
    cancers = list(meta.index[meta.disease_stage == 'Cancer'])
    controls = list(meta.index[meta.disease_stage != 'Cancer'])
    folds, selected = [], []
    for heldout in cancers:
        training = [s for s in cancers if s != heldout]
        if len(training) < 2 or len(controls) < 2:
            continue
        # All selection steps exclude held-out sample; expression eligibility is training-only.
        chosen, effect = select_cancer_panel(counts, logcpm, training, controls)
        if not chosen:
            continue
        score = logcpm[chosen].mean(axis=1)
        control_mean = float(score.loc[controls].mean())
        for gene in chosen:
            selected.append(dict(heldout_sample=heldout, gene=gene, training_logcpm_difference=float(effect[gene])))
        folds.append(dict(heldout_sample=heldout, signature_genes=';'.join(chosen),
                          heldout_score=float(score.loc[heldout]), training_cancer_mean=float(score.loc[training].mean()),
                          training_control_mean=control_mean, heldout_minus_control=float(score.loc[heldout] - control_mean)))
    pd.DataFrame(folds).to_csv(ROOT / 'results/signatures/lopo_cancer.csv', index=False)
    pd.DataFrame(selected).to_csv(ROOT / 'results/signatures/lopo_selected_genes.csv', index=False)
    if folds:
        f = pd.DataFrame(folds)
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.bar(f.heldout_sample, f.heldout_minus_control, color='#527ca5')
        ax.axhline(0, color='black'); ax.set_ylabel('Held-out score minus training control mean')
        ax.set_title('Cancer signature transfer: genes selected without held-out sample')
        save('progression/lopo_cancer')
    # Descriptive contrasts with a sample as replicate, no cell-level disease p-values.
    effects = pd.DataFrame(index=genes)
    for stage in STAGES[1:]:
        ids = meta.index[meta.disease_stage == stage]
        ref = meta.index[meta.disease_stage == STAGES[0]]
        effects[stage + '_minus_normal_logcpm'] = logcpm.loc[ids].mean() - logcpm.loc[ref].mean()
    effects.to_csv(ROOT / 'results/markers/epithelial_sample_level_effects.csv.gz')

def select_cancer_panel(counts, logcpm, training, controls):
    """Select only using explicitly supplied training rows, never held-out rows."""
    if set(training) & set(controls):
        raise ValueError('Training cancer and control samples overlap')
    eligible = (counts.loc[training + controls] >= 10).sum(axis=0) >= 3
    effect = (logcpm.loc[training].mean() - logcpm.loc[controls].mean())[eligible]
    return list(effect[effect > 0].nlargest(20).index), effect

def trajectory(e):
    from scipy.sparse.csgraph import connected_components
    labels = connected_components(e.obsp['connectivities'], directed=False)[1]
    counts = pd.crosstab(labels, e.obs.disease_stage)
    counts.to_csv(ROOT / 'results/trajectory_component_stages.csv')
    normal = np.flatnonzero(e.obs.disease_stage.to_numpy() == STAGES[0])
    if not len(normal):
        (ROOT / 'results/trajectory_status.json').write_text(json.dumps({'status': 'skipped', 'reason': 'No HPV-negative normal epithelial cells'}))
        return
    # Largest connected component that contains a normal cell; do not invent unreachable times.
    eligible_components = np.unique(labels[normal])
    sizes = np.bincount(labels)
    component = eligible_components[np.argmax(sizes[eligible_components])]
    t = e[labels == component].copy()
    t.obs['leiden'] = t.obs['leiden'].cat.remove_unused_categories()
    root_candidates = np.flatnonzero(t.obs.disease_stage.to_numpy() == STAGES[0])
    center = t.obsm['X_pca_harmony'][root_candidates].mean(axis=0)
    root = root_candidates[np.argmin(((t.obsm['X_pca_harmony'][root_candidates] - center) ** 2).sum(axis=1))]
    t.uns['iroot'] = int(root)
    sc.tl.diffmap(t, random_state=SEED)
    sc.tl.dpt(t)
    sc.tl.paga(t, groups='leiden')
    groups = t.obs.leiden.cat.categories
    pd.DataFrame(t.uns['paga']['connectivities'].toarray(), index=groups, columns=groups).to_csv(ROOT / 'results/epithelial_paga_connectivities.csv')
    sc.pl.paga(t, show=False, title='Epithelial PAGA: transcriptional connectivity, not lineage')
    save('progression/epithelial_paga')
    t.obs[['sample_id', 'disease_stage', 'leiden', 'dpt_pseudotime']].to_csv(ROOT / 'results/epithelial_dpt.csv.gz')
    embeddings(t, 'epithelial_trajectory', ['dpt_pseudotime'])
    status = dict(status='exploratory', root_cell=t.obs_names[root], included_cells=t.n_obs,
                  excluded_disconnected_cells=e.n_obs - t.n_obs, components=int(labels.max() + 1),
                  warning='Normal-rooted graph distance, not evidence of longitudinal progression; root and integration dependent')
    (ROOT / 'results/trajectory_status.json').write_text(json.dumps(status, indent=2))
