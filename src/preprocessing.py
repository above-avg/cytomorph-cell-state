import gc
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from .config import ROOT, SEED, QC
from .plotting import qc_plot, save

def load_qc(samples):
    retained, stats, metrics, sensitivity = [], [], [], []
    for row in samples.to_dict('records'):
        print('Loading', row['sample_id'], flush=True)
        a = sc.read_10x_mtx(ROOT / 'data/raw', prefix=row['prefix'] + '.', var_names='gene_symbols', cache=False)
        a.var_names_make_unique()
        if a.X.data.size == 0 or (a.X.data < 0).any() or not np.allclose(a.X.data, np.round(a.X.data)):
            raise ValueError('Expected nonnegative integer counts')
        for key in ['sample_id', 'patient_id', 'disease_stage', 'HPV_status', 'filename_histology']:
            a.obs[key] = row[key]
        a.obs['barcode'] = a.obs_names
        a.obs_names = row['sample_id'] + ':' + a.obs_names
        a.var['mt'] = a.var_names.str.startswith('MT-')
        if a.var.mt.sum() == 0:
            raise ValueError('No mitochondrial genes: gene symbols must be checked')
        sc.pp.calculate_qc_metrics(a, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
        keep = ((a.obs.n_genes_by_counts >= QC['min_genes']) & (a.obs.n_genes_by_counts <= QC['max_genes']) & (a.obs.pct_counts_mt <= QC['max_mt']))
        metrics.append(a.obs[['sample_id', 'n_genes_by_counts', 'total_counts', 'pct_counts_mt']].copy())
        stats.append(dict(sample_id=row['sample_id'], disease_stage=row['disease_stage'], input_cells=a.n_obs,
                          retained_cells=int(keep.sum()), retention_fraction=float(keep.mean()),
                          median_genes=float(a.obs.n_genes_by_counts.median()), median_counts=float(a.obs.total_counts.median()),
                          median_mt=float(a.obs.pct_counts_mt.median())))
        for max_genes, max_mt in [(5000, 10), (5000, 15), (5000, 20), (7500, 15)]:
            k = (a.obs.n_genes_by_counts >= 200) & (a.obs.n_genes_by_counts <= max_genes) & (a.obs.pct_counts_mt <= max_mt)
            sensitivity.append(dict(sample_id=row['sample_id'], max_genes=max_genes, max_mt=max_mt, retained_cells=int(k.sum())))
        a = a[keep].copy()
        if a.n_obs < 50:
            raise ValueError('Too few retained cells; inspect sample QC')
        a.write_h5ad(ROOT / f'data/processed/{row["sample_id"]}_counts.h5ad', compression='gzip')
        retained.append(a)
    pd.DataFrame(stats).to_csv(ROOT / 'results/qc_summary.csv', index=False)
    pd.DataFrame(sensitivity).to_csv(ROOT / 'results/qc_threshold_sensitivity.csv', index=False)
    qc_plot(pd.concat(metrics))
    merged = ad.concat(retained, join='inner', merge='same')
    del retained, metrics
    gc.collect()
    sc.pp.filter_genes(merged, min_cells=3)
    merged.write_h5ad(ROOT / 'data/processed/counts.h5ad', compression='gzip')
    return merged

def normalize(a):
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    a.raw = a
    sc.pp.highly_variable_genes(a, flavor='seurat', n_top_genes=2000, batch_key='sample_id')
    a.var.to_csv(ROOT / 'results/gene_features.csv')
    return a

def graph(a, prefix, resolution=0.6):
    """Centered sparse PCA; compare uncorrected graph with Harmony by sample."""
    import harmonypy as hm
    print('PCA', prefix, a.shape, flush=True)
    sc.pp.pca(a, n_comps=30, mask_var='highly_variable', svd_solver='arpack', random_state=SEED)
    sc.pl.pca_variance_ratio(a, log=True, show=False); save(f'umap/{prefix}_pca_variance')
    sc.pp.neighbors(a, n_neighbors=15, n_pcs=30, random_state=SEED)
    sc.tl.umap(a, random_state=SEED)
    a.obsm['X_umap_uncorrected'] = a.obsm['X_umap'].copy()
    # leidenalg avoids Scanpy 1.10's igraph RNG randint overflow on Windows.
    sc.tl.leiden(a, resolution=resolution, key_added='leiden_uncorrected', random_state=SEED, flavor='leidenalg', n_iterations=2, directed=False)
    a.obsp['uncorrected_connectivities'] = a.obsp['connectivities'].copy()
    print('Harmony', prefix, flush=True)
    # The package defaults derive 100 Harmony clusters for this cohort. A
    # smaller, explicitly fixed cluster count keeps the fit tractable while
    # retaining enough centroids to model sample-specific structure.
    harmony = hm.run_harmony(a.obsm['X_pca'], a.obs, 'sample_id', random_state=SEED,
                             nclust=50, max_iter_harmony=10, max_iter_kmeans=10,
                             verbose=False)
    a.obsm['X_pca_harmony'] = harmony.Z_corr.T.astype('float32')
    sc.pp.neighbors(a, use_rep='X_pca_harmony', n_neighbors=15, random_state=SEED)
    sc.tl.umap(a, random_state=SEED)
    sc.tl.leiden(a, resolution=resolution, key_added='leiden', random_state=SEED, flavor='leidenalg', n_iterations=2, directed=False)
    return a
