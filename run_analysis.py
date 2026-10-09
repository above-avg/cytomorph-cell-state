"""Run the real cohort end to end; checkpoints support inspection and resumption."""
import argparse
import gc
import json
import platform
import time
from importlib.metadata import version
import pandas as pd
import scanpy as sc
from src.config import ROOT, SEED, QC
from src.download import prepare
from src.preprocessing import load_qc, normalize, graph
from src.analysis import annotate, markers, epithelial, immune, integration_diagnostics, pseudobulk_and_lopo, trajectory

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--resume', action='store_true', help='Use existing checkpoints; remove stale ones after changing parameters')
    parser.add_argument('--stage', choices=['all', 'qc', 'cluster', 'epithelial', 'robustness'], default='all')
    args = parser.parse_args()
    start = time.time()
    for directory in ['data/processed', 'results/markers', 'results/signatures', 'figures']:
        (ROOT / directory).mkdir(parents=True, exist_ok=True)
    # One job avoids multiprocessing pipe restrictions in Windows sandboxes.
    sc.settings.n_jobs = 1
    sc.settings.verbosity = 2
    if args.resume and (ROOT / 'data/metadata/samples.csv').exists():
        samples = pd.read_csv(ROOT / 'data/metadata/samples.csv')
    else:
        samples = prepare()
    final = ROOT / 'data/processed/annotated.h5ad'
    epi = ROOT / 'data/processed/epithelial.h5ad'
    if args.stage in ['all', 'qc', 'cluster']:
        if args.resume and final.exists():
            a = sc.read_h5ad(final)
        else:
            countfile = ROOT / 'data/processed/counts.h5ad'
            a = sc.read_h5ad(countfile) if args.resume and countfile.exists() else load_qc(samples)
            if args.stage == 'qc':
                return
            normalize(a); graph(a, 'all_cells'); annotate(a)
            integration_diagnostics(a, 'all_cells')
            markers(a, 'leiden', 'all_cells')
            a.write_h5ad(final, compression='gzip')
        if args.stage == 'cluster':
            return
    elif args.stage == 'epithelial':
        a = sc.read_h5ad(final)
    if args.stage in ['all', 'epithelial']:
        immune(a)
        e = sc.read_h5ad(epi) if args.resume and epi.exists() else epithelial(a)
        del a
        gc.collect()
    elif args.stage == 'robustness':
        e = sc.read_h5ad(epi)
    else:
        return
    pseudobulk_and_lopo(e, samples)
    trajectory(e)
    manifest = dict(dataset='GSE208653', seed=SEED, qc=QC, resumed=args.resume, python=platform.python_version(),
                    versions={p: version(p) for p in ['scanpy', 'anndata', 'numpy', 'scipy', 'pandas', 'igraph', 'leidenalg', 'harmonypy']},
                    elapsed_seconds=round(time.time() - start, 1), epithelial_cells=e.n_obs,
                    samples=samples.sample_id.tolist(), status='completed')
    (ROOT / 'results/run_manifest.json').write_text(json.dumps(manifest, indent=2))
    print(json.dumps(manifest, indent=2), flush=True)

if __name__ == '__main__':
    main()
