# Analysis decisions

## Question and design

Do epithelial states and immune composition differ across HPV-negative normal,
HPV-positive normal, HSIL, and cancer tissue? Which patterns recur across samples?
These groups form a cross-sectional clinical ordering, not a longitudinal record
of the same patients. Disease group is confounded with HPV status for several
contrasts. The two normal groups provide the direct HPV-status comparison.

## QC

Each supplied sample is loaded independently. Gene symbols are made unique;
barcodes receive sample prefixes. Counts must be nonnegative integers, and `MT-`
genes must exist. All per-cell metadata remain attached through subsetting.

Primary thresholds are prespecified exploratory rules: 200–5,000 detected genes
and mitochondrial counts at most 15%. The workflow exports pre-filter distributions,
counts-versus-genes plots, retention per sample, and alternative retention for 10%
and 20% mitochondrial thresholds and a 7,500-gene ceiling. These thresholds are
not claimed to have been optimized after looking at the result. Sensitivity tables
measure **retention only**, not full downstream stability under every threshold.
The gene ceiling is not a substitute for doublet detection. Doublets and ambient
RNA remain unmodeled limitations; the data are already supplied count matrices.

## Normalization and representation

Retained raw counts are saved per sample before any normalization. The pooled
matrix uses shared genes detected in at least three cells. Counts are normalized
to 10,000 per cell and log1p transformed. Full normalized expression is preserved
in `.raw` for markers. Seurat-flavor highly variable genes are selected per batch,
combining evidence across samples, with a target of 2,000 genes.

Thirty centered sparse principal components are computed on HVGs without unit
variance scaling or regression. This choice avoids densifying the full matrix;
high-variance genes can have more influence. Cell cycle is retained as possible
biology. A 15-neighbor graph, UMAP, and Leiden clustering are first computed without
correction. Harmony then adjusts PCA coordinates by sample and the graph is rebuilt.
The fit uses 50 centroids, at most 10 Harmony rounds, and at most 10 k-means
rounds to keep the nine-sample analysis tractable on a desktop machine.
Seed: 42. Leiden uses the `leidenalg` backend with two iterations, an undirected
graph, and resolution 0.6 for all cells and 0.5 for epithelial cells. This backend
avoids the igraph RNG integer overflow observed with Scanpy 1.10 on Windows.

Harmony can remove disease-associated variation because each sample belongs to
one disease group. Both embeddings and cluster-by-sample tables are exported.
Adjusted Rand index compares baseline and corrected cluster assignments; it is
not an integration quality score. A low dominant-sample fraction alone does not
prove successful integration. Expression-based analyses use uncorrected counts
or log-normalized expression, not Harmony coordinates.

## Broad annotation and subset analysis

Canonical marker panels are declared in `src/config.py`. Scanpy control-matched
gene scores are averaged by corrected Leiden cluster. The highest scoring panel
assigns a provisional label only if its score is positive and its margin over
the next panel is at least 0.1; otherwise the label is Ambiguous. This heuristic
is **not reference-based or expert-validated annotation**. The marker dotplot,
score table, competing label, and cluster markers are exported for review.
T and NK are intentionally combined; smooth-muscle/pericyte resolution is limited.

Epithelial cells are reclustered using new subset-specific HVGs, PCA, Harmony,
neighbors, UMAP, and Leiden. A cluster number denotes a transcriptional grouping,
not malignancy or disease stage. Composition is exported by sample and by stage;
stage summaries average sample fractions equally instead of pooling all cells.
These are recovered-cell fractions, not estimates of in vivo tissue abundance.

## Marker rankings, signatures, and replication

Cluster-versus-rest Welch t-test rankings are descriptive cell-level summaries over
the 500 most variable genes among each stage's 2,000 selected highly variable genes, with Scanpy's adjusted p-values retained for transparency. They must not be interpreted
as patient-level disease evidence. Selecting and testing clusters on the same
expression data also makes these p-values exploratory.

The published epithelial marker pairs (SLC5A8/DERL3, VSIG1/CASC9, CASP14/CALML5)
are hypotheses to inspect, not labels assigned to the new clusters. Prespecified
small panels for cell cycle, keratinization, interferon response, cytotoxicity,
and inhibitory receptors are also inspected. Scores are arithmetic means of
log-normalized expression over present genes. They are not enrichment tests or
validated biomarkers. Coverage is exported, including missing genes. Receptor
expression alone does not establish immune exhaustion. T/NK-specific summaries
help distinguish changing immune composition from changing expression within T/NK.

Raw epithelial counts are summed by sample to form pseudobulk counts. Log1p CPM
is used for descriptive stage differences. No formal patient-level differential
expression p-values are claimed in this small cohort.

For each held-out cancer sample, a 20-gene cancer panel is selected from **only
the other cancer samples and all noncancer samples** using training mean logCPM
differences. Training-only eligibility requires at least 10 counts in at least
three training samples. The panel is scored on the held-out sample and compared
with the training noncancer mean. Gene selection excludes the held-out sample,
but upstream annotation/clustering used the full cohort. Thus this is conditional
sample-level transfer of an expression panel, **not fully nested validation of
the entire pipeline**, independent validation, or a predictive classifier.

## Exploratory trajectory

DPT uses the largest Harmony graph component containing HPV-negative normal
epithelial cells. Its root is the normal cell closest to the normal PCA centroid.
Disconnected cells are excluded and counted, not assigned fabricated pseudotime.
PAGA connectivity is computed for subset clusters. DPT is normal-rooted graph
distance; it cannot prove normal-to-cancer transformation. Tissue type, histology,
sample effects, root choice, and integration can all shape the result.

## Scope

No paired PAP-smear analysis, inferred CNV malignancy calling, clinical diagnostic
validation, causal HPV attribution, or cell-cell communication analysis is claimed.
This is an exploratory reanalysis of a small public cohort, supported by actual
executed outputs rather than fabricated findings.
