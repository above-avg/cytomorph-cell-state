# Questions and execution notes

This log records the questions, implementation decisions, and unresolved issues
that shaped the analysis. It distinguishes results from interpretations and notes
where additional biological review or data would be useful.

| Question or issue | What was done | What remains unresolved |
|---|---|---|
| Can the PBMC tutorial answer the cervical question? | Replaced the initial PBMC idea with the requested GSE208653 cohort before running analysis. | Tissue-specific annotation needs stronger reference support. |
| Is the sample mapping trustworthy? | Parsed disease state and filenames from the complete GEO SOFT record instead of copying a proposed table. | Separate donor IDs and detailed clinical covariates are absent from the exported mapping. |
| What does “raw” mean here? | Confirmed the archive contains sparse count matrices with barcodes and gene features. | Alignment, ambient RNA and original cell calling cannot be evaluated as if FASTQ had been processed here. |
| Why did initial QC fail? | `calculate_qc_metrics` returns results unless `inplace=True`; corrected the call and reran before any results were accepted. | API defaults need explicit attention when translating tutorials. |
| Why did graph/clustering execution fail on Windows? | Local worker pipes were blocked by the sandbox. Pinned numerical versions, then isolated an integer overflow to Scanpy's igraph random-state adapter during Leiden. Switched to the supported `leidenalg` backend and resumed from saved counts. | Cross-platform numerical outputs can differ slightly even with a fixed seed. |
| Are high gene counts definitely doublets? | Kept a transparent gene ceiling and reported retention sensitivity. | A dedicated doublet/ambient-RNA analysis is still needed. |
| Does a well-mixed UMAP prove integration worked? | Saved baseline and Harmony embeddings plus sample contributions to each cluster. | Sample correction may suppress stage biology; visual mixing cannot resolve this by itself. |
| Are canonical marker panels enough to name every cluster? | Exported competing scores, margins, marker dotplots and an Ambiguous option. | Labels are provisional; expert review or reference mapping remains valuable. |
| Do many cells imply many biological replicates? | Used sample-level fractions, expression summaries, pseudobulk and held-out cancer panels. | Only 2–3 samples per group; uncertainty remains large. |
| Can published marker pairs be treated as confirmed findings? | Inspected them as prespecified hypotheses and exported gene coverage. | Patterns may differ by patient, histology, annotation or processing. |
| Does pseudotime establish cancer progression? | Restricted DPT to a connected component with a declared normal root. | Cross-sectional graph distance is not longitudinal lineage or causality. |

The main biological observation and open question are written in the README
only after examining the completed result tables.
