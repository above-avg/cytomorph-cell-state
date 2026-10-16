"""Create notebook walkthroughs that inspect the saved full analysis."""
import sys
from pathlib import Path
import nbformat as nbf

ROOT = Path(__file__).resolve().parent
BOOT = '''from pathlib import Path
import sys, json
import pandas as pd
from IPython.display import display, Image
ROOT = Path.cwd()
if not (ROOT / 'run_analysis.py').exists():
    ROOT = ROOT.parent
sys.path.insert(0, str(ROOT))
def table(relative, **kwargs):
    return pd.read_csv(ROOT / relative, **kwargs)
def figure(relative):
    display(Image(filename=str(ROOT / 'figures' / relative)))
'''

def md(s): return nbf.v4.new_markdown_cell(s)
def code(s): return nbf.v4.new_code_cell(s)

SPECS = [
    ('01_data_loading', 'The cohort and its biological question', [
        md('Do epithelial states differ across cervical disease groups, and are patterns shared across samples? The clinical groups are cross-sectional. No paired PAP-smear images are available.\n\nThe full analysis is performed by `run_analysis.py`; these walkthrough notebooks inspect its saved real outputs. Download and checkpoint commands are included for reproduction, without repeating expensive processing in every notebook.'),
        code("samples = table('data/metadata/samples.csv')\ndisplay(samples[['sample_id', 'disease_stage', 'HPV_status', 'filename_histology']])\nprint(samples.groupby('disease_stage').size())"),
        md('Sample accessions act as replicate proxies. Separate donor identifiers are not supplied by GEO. Cancer filenames distinguish SCC and ADC; this can complicate a shared cancer signature.'),
        code("import scanpy as sc\na = sc.read_h5ad(ROOT / 'data/processed/GSM6360680_counts.h5ad', backed='r')\nprint(a)\ndisplay(a.obs[['barcode', 'sample_id', 'patient_id', 'disease_stage', 'HPV_status']].head())\na.file.close()"),
        md('Reproduce data acquisition with `python -m src.download`. The download code validates archive paths and writes checksums. Count matrices are deposited expression data, not FASTQ.\n\n**Question:** What additional clinical covariates would help separate disease, histology, and donor effects?'),
    ]),
    ('02_quality_control', 'QC per sample, before pooling', [
        md('The filters are prespecified exploratory thresholds: 200–5,000 genes and ≤15% mitochondrial counts. The pre-filter distributions are inspected below. Do not rewrite this as retrospective threshold optimization.'),
        code("qc = table('results/qc_summary.csv')\ndisplay(qc)\nprint('Input:', qc.input_cells.sum(), 'Retained:', qc.retained_cells.sum())"),
        code("figure('qc/per_sample.png')\nfigure('qc/counts_vs_genes.png')"),
        code("sensitivity = table('results/qc_threshold_sensitivity.csv')\ndisplay(sensitivity.pivot(index='sample_id', columns=['max_genes', 'max_mt'], values='retained_cells'))"),
        md('The alternative thresholds evaluate retention, not full reclustering stability. High gene count filtering does not establish doublet removal.\n\nRecompute this stage with `python run_analysis.py --stage qc`. Counts are preserved in per-sample checkpoints before normalization.\n\n**Question:** Why might a fixed mitochondrial threshold remove different proportions from different tissues?'),
    ]),
    ('03_integration_clustering', 'Normalization, graph construction, and sample correction', [
        md('The CLI normalizes to 10,000 counts per cell, applies log1p, selects 2,000 batch-aware HVGs, and computes 30 centered sparse PCs. A 15-neighbor graph supports seeded UMAP and Leiden. Harmony adjusts sample effects, followed by a new graph. Both results are retained.\n\n`python run_analysis.py --resume --stage cluster` runs this stage from saved counts. Delete stale checkpoints when changing parameters.'),
        code("from src.preprocessing import normalize, graph\nimport inspect\nprint(inspect.getsource(normalize))\nprint(inspect.getsource(graph))"),
        code("display(json.loads((ROOT / 'results/all_cells_integration_diagnostics.json').read_text()))\nfigure('umap/all_cells_baseline_sample_id.png')\nfigure('umap/all_cells_sample_id.png')"),
        code("figure('umap/all_cells_baseline_disease_stage.png')\nfigure('umap/all_cells_disease_stage.png')"),
        md('Sample mixing is not proof that biology is preserved. Stage and sample are linked, so Harmony may overcorrect. The adjusted Rand index describes agreement between cluster assignments, not correctness. Expression summaries below always use counts or uncorrected log-normalized expression.\n\n**Question:** Which apparent stage separation changes most after correction, and is that change technical or biological?'),
    ]),
    ('04_cell_type_annotation', 'Provisional broad identities from marker evidence', [
        code("from src.config import MARKERS\ndisplay(MARKERS)\ndisplay(table('results/annotation_review.csv'))"),
        code("figure('annotation/broad_marker_dotplot.png')\nfigure('umap/all_cells_cell_type.png')"),
        md('Control-matched marker-panel scores are averaged by cluster. The best positive score assigns a label only with margin ≥0.1. Ambiguous clusters are retained. This heuristic needs marker review and is not expert-validated annotation. T/NK is a combined category; tumor malignancy is not inferred from epithelial expression.'),
        code("display(table('results/markers/all_cells_top20.csv').groupby('group').head(5))\nfigure('progression/broad_composition_by_sample.png')"),
        md('Cell-level cluster marker p-values are descriptive after clustering on the same expression matrix; they are not evidence that disease effects are replicated across patients.\n\n**Question:** Do the competing label and top markers support each provisional assignment?'),
    ]),
    ('05_epithelial_progression', 'Epithelial states and the cytology connection', [
        md('Epithelial cells are selected from provisional broad labels, then given new HVGs, PCA, graphs, UMAP and Leiden clusters. Cluster numbers are not stage names.\n\nCytology shows cellular morphology; transcriptomics measures molecular heterogeneity. These samples are unpaired, so no cluster is mapped directly to a PAP-smear appearance.'),
        code("figure('umap/epithelial_leiden.png')\nfigure('umap/epithelial_disease_stage.png')\nfigure('umap/epithelial_sample_id.png')"),
        code("figure('progression/epithelial_composition_by_sample.png')\nfigure('progression/epithelial_composition_by_stage.png')"),
        code("display(table('results/epithelial_composition_sample_counts.csv'))\ndisplay(table('results/epithelial_composition_stage_mean_fractions.csv'))"),
        md('Stage bars average within-sample fractions equally; large samples do not dominate. Fractions concern recovered cells, not tissue abundance. A state enriched in one cancer sample is not automatically cancer-associated across the cohort.\n\n**Question:** Which epithelial state is shared by multiple cancer samples, and which is dominated by one sample?'),
    ]),
    ('06_differential_expression', 'Marker rankings, programs, and sample-level summaries', [
        md('Cluster markers use Welch t-test rankings versus the remaining cells, restricted to the 500 most variable genes among the selected HVGs. For disease-group contrasts, epithelial raw counts are summed by sample into pseudobulk. LogCPM differences are descriptive effect summaries; no patient-level significance test is claimed.'),
        code("display(table('results/markers/epithelial_top20.csv').groupby('group').head(5))\ndisplay(table('results/signatures/epithelial_pseudobulk_samples.csv'))"),
        code("coverage = table('results/signatures/epithelial_gene_coverage.csv')\ndisplay(coverage)\nfigure('progression/epithelial_signature_scores.png')"),
        code("effects = table('results/markers/epithelial_sample_level_effects.csv.gz', index_col=0)\ndisplay(effects.sort_values('Cancer_minus_normal_logcpm', ascending=False).head(20))"),
        code("figure('progression/immune_composition_by_sample.png')\nfigure('progression/T_NK_signature_scores.png')"),
        md('Published epithelial marker pairs are hypotheses, not confirmed labels. Small gene-panel scores are mean log-normalized expression, not enrichment tests. Inhibitory-receptor expression cannot alone prove exhaustion. Compare T/NK-specific expression with immune composition.\n\n**Question:** Could a stage-level signature be explained by different cell composition or cancer histology?'),
    ]),
    ('07_patient_robustness', 'Held-out samples, uncertainty, and exploratory trajectory', [
        md('Each cancer sample is withheld while a 20-gene panel is selected using the other two cancer samples and all noncancer samples. Eligibility and selection use training samples only. Upstream annotation used the full cohort; this is conditional panel transfer, not fully nested validation or an independent test cohort.'),
        code("folds = table('results/signatures/lopo_cancer.csv')\ndisplay(folds)\nfigure('progression/lopo_cancer.png')"),
        code("selected = table('results/signatures/lopo_selected_genes.csv')\ndisplay(selected.groupby('gene').heldout_sample.nunique().sort_values(ascending=False).head(20))"),
        md('A positive held-out contrast suggests transfer within this small cohort. A failed fold is evidence of heterogeneity, not a reason to discard that sample. The unit test changes held-out expression drastically and verifies that the selected training panel is unchanged.'),
        code("display(table('results/signatures/lopo_selected_genes.csv').groupby('heldout_sample').gene.apply(list))"),
        code("display(json.loads((ROOT / 'results/trajectory_status.json').read_text()))\nfigure('umap/epithelial_trajectory_dpt_pseudotime.png')"),
        md('Normal-rooted DPT is graph distance on one connected component. It does not show a cell becoming cancerous. Root, graph, sample and histology dependence remain unresolved.\n\n**Next experiment:** repeat annotation and QC sensitivity, then validate a prespecified expression panel in an independent cohort with more biological replicates.'),
    ]),
]

def main():
    out = ROOT / 'notebooks'
    out.mkdir(exist_ok=True)
    completed = []
    for name, title, cells in SPECS:
        nb = nbf.v4.new_notebook(cells=[md('# ' + title), code(BOOT), *cells])
        nb.metadata.kernelspec = dict(name='python3', display_name='Python 3', language='python')
        nb.metadata.language_info = dict(name='python', version='3.10')
        print('Writing', name, flush=True)
        nbf.validate(nb)
        nbf.write(nb, out / (name + '.ipynb'))
        completed.append(name)
    import json
    (ROOT / 'results/notebook_manifest.json').write_text(json.dumps({'status': 'created', 'notebooks': completed}, indent=2))

if __name__ == '__main__':
    main()
