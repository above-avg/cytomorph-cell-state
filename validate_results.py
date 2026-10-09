"""Check biological metadata, accounting, and generated analysis artifacts."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
from src.config import ROOT, STAGES, QC

def main():
    samples = pd.read_csv(ROOT / 'data/metadata/samples.csv')
    assert len(samples) == 9 and samples.sample_id.is_unique
    assert samples.disease_stage.value_counts().to_dict() == dict(zip(STAGES, [2, 2, 2, 3]))
    qc = pd.read_csv(ROOT / 'results/qc_summary.csv')
    assert len(qc) == 9 and (qc.retained_cells <= qc.input_cells).all()
    a = sc.read_h5ad(ROOT / 'data/processed/annotated.h5ad', backed='r')
    assert a.n_obs == qc.retained_cells.sum() and a.obs_names.is_unique
    assert set(a.obs.sample_id) == set(samples.sample_id)
    assert np.isfinite(a.obsm['X_umap']).all()
    assert (a.obs.n_genes_by_counts.between(QC['min_genes'], QC['max_genes'])).all()
    assert (a.obs.pct_counts_mt <= QC['max_mt']).all()
    for key in ['disease_stage', 'patient_id', 'HPV_status', 'cell_type']:
        assert a.obs[key].notna().all()
    a.file.close()
    e = sc.read_h5ad(ROOT / 'data/processed/epithelial.h5ad', backed='r')
    assert (e.obs.cell_type == 'Epithelial').all() and e.obs_names.is_unique
    assert np.isfinite(e.obsm['X_umap']).all()
    pbmeta = pd.read_csv(ROOT / 'results/signatures/epithelial_pseudobulk_samples.csv')
    assert pbmeta.cells.sum() == e.n_obs
    e.file.close()
    for name in ['broad_composition', 'epithelial_composition', 'immune_composition']:
        f = pd.read_csv(ROOT / f'results/{name}_sample_fractions.csv', index_col=0)
        assert np.allclose(f.sum(axis=1), 1) and np.isfinite(f.to_numpy()).all()
    folds = pd.read_csv(ROOT / 'results/signatures/lopo_cancer.csv')
    assert len(folds) == 3 and folds.heldout_sample.is_unique
    assert np.isfinite(folds.select_dtypes('number')).all().all()
    count = pd.read_csv(ROOT / 'results/signatures/epithelial_pseudobulk_counts.csv.gz', index_col=0)
    assert (count.to_numpy() >= 0).all() and np.allclose(count, np.round(count))
    assert json.loads((ROOT / 'results/run_manifest.json').read_text())['status'] == 'completed'
    pngs = list((ROOT / 'figures').rglob('*.png'))
    assert len(pngs) >= 20
    report = dict(status='passed', samples=len(samples), input_cells=int(qc.input_cells.sum()),
                  retained_cells=int(qc.retained_cells.sum()), epithelial_cells=int(pbmeta.cells.sum()),
                  lopo_folds=len(folds), generated_figures=len(pngs))
    (ROOT / 'results/validation.json').write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))

if __name__ == '__main__':
    main()
