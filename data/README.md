# Data provenance

The analysis uses **GSE208653**, the public cervical-tissue single-cell count matrices
from Guo and colleagues. These are deposited expression matrices, not raw FASTQ
reads. This repository does not redo alignment, tissue sampling, or cell calling.

- [GEO series](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE208653)
- [Original study](https://pmc.ncbi.nlm.nih.gov/articles/PMC10040725/)
- [Expression archive](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE208nnn/GSE208653/suppl/GSE208653_RAW.tar)
- [Full GEO SOFT metadata](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE208nnn/GSE208653/soft/GSE208653_family.soft.gz)

Run `python -m src.download` from the repository root. This downloads the archive
(643,706,880 bytes), validates flat archive paths, extracts 27 matrix/barcode/feature
files, builds `metadata/samples.csv` from GEO metadata, and records SHA-256 hashes.
Allow several GB of disk space for the environment and analysis checkpoints.

The metadata distinguishes 2 HPV-negative normal, 2 HPV-positive normal, 2 HSIL,
and 3 cancer samples. `patient_id` uses the GSM accession as a **sample replicate
proxy**, because GEO does not provide a separate donor identifier. The cancer
filenames encode SCC for GSM6360686/687 and ADC for GSM6360688; this is retained as
`filename_histology` rather than asserted as a newly validated pathology diagnosis.

`raw/` holds downloaded matrices and `processed/` holds local AnnData checkpoints.
Both are ignored by Git. `metadata/` contains the small source metadata, accession
mapping, and checksums. No PAP images or paired cytology/transcriptome data are used.
