"""Download public GEO matrices without putting large files into Git."""
import gzip
import hashlib
import json
import re
import tarfile
import urllib.request
from pathlib import Path
import pandas as pd
from .config import ROOT, STAGES

BASE = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE208nnn/GSE208653'

def download(url, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.exists():
        tmp = path.with_suffix(path.suffix + '.part')
        urllib.request.urlretrieve(url, tmp)
        tmp.replace(path)

def prepare():
    meta = ROOT / 'data/metadata'
    soft = meta / 'GSE208653_family.soft.gz'
    download(f'{BASE}/soft/GSE208653_family.soft.gz', soft)
    text = gzip.open(soft, 'rt').read()
    rows = []
    for block in text.split('^SAMPLE = ')[1:]:
        accession = block.splitlines()[0].strip()
        title = re.search(r'!Sample_title = (.*)', block).group(1)
        state = re.search(r'disease state: (.*)', block).group(1).strip()
        files = re.findall(r'!Sample_supplementary_file_\d+ = (.*)', block)
        prefixes = {Path(f.strip()).name.split('.barcodes')[0].split('.features')[0].split('.matrix')[0] for f in files}
        assert len(prefixes) == 1
        prefix = prefixes.pop()
        stage = dict(NO_HPV=STAGES[0], N_HPV=STAGES[1], HSIL_HPV=STAGES[2], CA_HPV=STAGES[3])[state]
        rows.append(dict(sample_id=accession, patient_id=accession,
                         patient_id_note='sample accession used as replicate proxy; no donor ID in GEO',
                         disease_stage=stage, HPV_status='Negative' if state == 'NO_HPV' else 'Positive',
                         sample_title=title, prefix=prefix,
                         filename_histology='ADC' if '_ADC_' in prefix else ('SCC' if '_SCC_' in prefix else 'not_applicable')))
    frame = pd.DataFrame(rows)
    assert len(frame) == 9 and frame.sample_id.is_unique
    frame.to_csv(meta / 'samples.csv', index=False)
    archive = ROOT / 'data/raw/GSE208653_RAW.tar'
    download(f'{BASE}/suppl/GSE208653_RAW.tar', archive)
    # Never trust archive paths: only extract ordinary flat matrix files.
    with tarfile.open(archive) as tar:
        for member in tar.getmembers():
            if not member.isfile() or Path(member.name).name != member.name:
                raise ValueError(f'Unexpected archive member: {member.name}')
            if not member.name.endswith(('.barcodes.tsv.gz', '.features.tsv.gz', '.matrix.mtx.gz')):
                raise ValueError(f'Unexpected matrix file: {member.name}')
            target = archive.parent / member.name
            if not target.exists():
                with tar.extractfile(member) as src, target.open('wb') as out:
                    import shutil
                    shutil.copyfileobj(src, out)
    hashes = {}
    for path in [soft, *sorted(archive.parent.glob('*.gz'))]:
        h = hashlib.sha256()
        with path.open('rb') as f:
            for chunk in iter(lambda: f.read(1024 * 1024), b''):
                h.update(chunk)
        hashes[path.relative_to(ROOT).as_posix()] = h.hexdigest()
    (meta / 'sha256.json').write_text(json.dumps(hashes, indent=2))
    return frame

if __name__ == '__main__':
    prepare()
