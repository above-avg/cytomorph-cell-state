from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
STAGES = ['Normal HPV-', 'Normal HPV+', 'HSIL', 'Cancer']
SEED = 42
# Prespecified exploratory thresholds; distributions and sensitivity are exported.
QC = dict(min_genes=200, max_genes=5000, max_mt=15)
MARKERS = {
    'Epithelial': ['EPCAM', 'KRT8', 'KRT18', 'KRT19', 'KRT5', 'KRT14'],
    'T/NK': ['CD3D', 'CD3E', 'TRAC', 'NKG7', 'GNLY'],
    'B': ['MS4A1', 'CD79A', 'CD79B', 'CD37'],
    'Plasma': ['MZB1', 'JCHAIN', 'SDC1', 'IGKC'],
    'Myeloid': ['LYZ', 'LST1', 'TYROBP', 'FCER1G', 'CTSS'],
    'Fibroblast': ['COL1A1', 'COL1A2', 'DCN', 'LUM'],
    'Endothelial': ['PECAM1', 'VWF', 'KDR', 'EMCN'],
    'Smooth muscle': ['ACTA2', 'TAGLN', 'MYH11', 'CNN1'],
    'Mast': ['TPSAB1', 'TPSB2', 'KIT', 'CPA3'],
}
SIGNATURES = {
    'HPV_normal_published': ['SLC5A8', 'DERL3'],
    'HSIL_published': ['VSIG1', 'CASC9'],
    'Cancer_published': ['CASP14', 'CALML5'],
    'Cell_cycle': ['MKI67', 'TOP2A', 'PCNA', 'MCM2', 'MCM5', 'TYMS'],
    'Keratinization': ['KRT1', 'KRT10', 'IVL', 'LOR', 'FLG', 'SPRR1B'],
    'Interferon_response': ['ISG15', 'IFIT1', 'IFIT3', 'MX1', 'OAS1', 'STAT1'],
    'Cytotoxicity': ['NKG7', 'GNLY', 'PRF1', 'GZMB', 'GZMH'],
    'Inhibitory_receptors': ['PDCD1', 'LAG3', 'TIGIT', 'HAVCR2', 'CTLA4'],
}
