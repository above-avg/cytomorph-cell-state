# Alignment with the cervical cancer multimodal project

Target: **Multimodal single-cell characterization of initiation and progression
of cervical cancer**, PI Saket Choudhary, as described in the supplied project brief.
This repository is preparatory work; it is not a claim of participation in the
professor's study or access to its private data.

## Working milestone

The implemented host-response extension compares HPV status within normal tissue
and subsequent lesion categories **within individual cell types**, using samples
as replicates. Its [demonstration report](host_response_results.md) provides actual
figures, sample-influence calculations, histology sensitivity, and a viral-feature
audit. It establishes a reproducible host-transcriptome baseline for later viral
and chromatin measurements.

| Project objective | Demonstrated now | Required next |
|---|---|---|
| Cell-type-specific HPV–host relationships | Host-program comparisons across tissue HPV-status groups; sample influence | Viral RNA measurements with reference and genotype provenance; infection cannot be assigned to a cell from tissue status |
| CIN initiation and tumour progression | Normal HPV-, normal HPV+, HSIL, cancer comparisons | Separate CIN grades and early/advanced cancer metadata; existing groups must not be relabelled as those grades |
| RNA and open-chromatin integration | Reproducible RNA baseline and preserved sample/cell identifiers | scATAC fragments or peak matrix, genome build, QC, and a validated pairing map or explicitly unpaired study design |
| Viral integration by cell type | Audit of whether current features explicitly represent HPV | Sequence-level host-virus junction evidence with cell barcode support and independent validation; expression alone does not locate integration |
| Indian-donor analysis | No claim of Indian-donor provenance | Authorised cohort access, donor provenance, clinical covariates, and sufficient biological replicates |
| Pan-cancer integration-site repository | No integration sites are fabricated from expression data | Curated evidence-backed site records across cancers, consistent coordinates, source provenance and validation status |

## Next bounded implementation

1. Inventory a permitted sequence-level dataset and its assay, reference build,
   HPV reference/genotype, barcode/UMI information, and clinical metadata.
2. Reprocess against documented host and viral references to distinguish host
   expression, viral expression, and candidate host-virus junction support.
3. Store candidate evidence with sample ID, barcode where supported, host build,
   host/viral contigs and breakpoint coordinates, coordinate convention, strands,
   read/UMI support, mapping-quality evidence, source accession and validation
   status. Retain negative and ambiguous results and guard against duplicate reads.
4. Benchmark candidate detection with known positive and negative controls before
   associating candidates with cell types. RNA-derived candidates require further
   validation before being presented as genomic integration events.
5. Add chromatin data only after confirming which cells or samples are paired.
   Report accessibility and expression associations separately from causal claims.

## Data boundaries

The current [GEO series](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE208653)
is the scRNA-seq series used in the existing analysis. Its description discusses
spatial transcriptomics; that does not establish paired open-chromatin data in
this checkout. The feature audit is a label search, not re-alignment to a viral
reference. No new private cohort, viral integration caller, scATAC analysis or
pan-cancer collection is represented as completed.
