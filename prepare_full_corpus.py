"""Resumable, parallel full-dataset SAFE conversion and NPE ring/atom census."""
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
import csv
import hashlib
import json
import multiprocessing as mp
from pathlib import Path
import time

PROJECT = Path(__file__).resolve().parent
OUTPUTS = PROJECT.parent
WORK = OUTPUTS.parent / 'work'
CACHE = WORK / 'chembl37_tokenizer_cache'
SAFE = None

def ring_systems(mol):
    rings = [set(ring) for ring in mol.GetRingInfo().AtomRings()]
    systems = []
    for ring in rings:
        connected = [old for old in systems if ring & old]
        combined = set(ring)
        for old in connected:
            combined.update(old)
            systems.remove(old)
        systems.append(combined)
    return systems

def prepare_batch(job):
    number, rows = job
    digest = hashlib.sha256(json.dumps(rows, separators=(',', ':')).encode()).hexdigest()
    prefix = CACHE / f'{number:05d}'
    stat_path = prefix.with_suffix('.json')
    if stat_path.exists():
        stats = json.loads(stat_path.read_text(encoding='utf-8'))
        if stats.get('input_sha256') == digest and stats.get('schema') == 1:
            return number, stats
    global SAFE
    if SAFE is None:
        from molecular_tokenizer.backends import SAFEBackend
        from rdkit import RDLogger
        RDLogger.DisableLog('rdApp.*')
        SAFE = SAFEBackend()
    from rdkit import Chem
    from molecular_tokenizer._vendor.graphbpe_utils import get_sub_molecule, molecule_to_smiles
    safe_count = npe_count = 0
    atoms, rings = Counter(), Counter()
    files = {name: prefix.with_suffix('.' + name) for name in ['safe.txt', 'npe.smi', 'rejects.jsonl']}
    temp = {name: Path(str(path) + '.tmp') for name, path in files.items()}
    with temp['safe.txt'].open('w', encoding='utf-8', newline='\n') as safe_out, temp['npe.smi'].open('w', encoding='utf-8', newline='\n') as npe_out, temp['rejects.jsonl'].open('w', encoding='utf-8', newline='\n') as rejected:
        for row_number, cid, smiles in rows:
            errors = {}
            try:
                safe_out.write(SAFE.representation(smiles) + '\n')
                safe_count += 1
            except Exception as error:
                errors['safe'] = f'{type(error).__name__}: {error}'
            try:
                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    raise ValueError('RDKit could not parse SMILES')
                Chem.Kekulize(mol, True)
                mol_rings = set()
                for system in ring_systems(mol):
                    ring = get_sub_molecule(mol, list(system), True)
                    mol_rings.add(molecule_to_smiles(ring, True))
                atoms.update(atom.GetSymbol() for atom in mol.GetAtoms())
                rings.update(mol_rings)  # presence per molecule, matching GraphBPE initialization
                npe_out.write(smiles + '\n')
                npe_count += 1
            except Exception as error:
                errors['npe'] = f'{type(error).__name__}: {error}'
            if errors:
                rejected.write(json.dumps({'row': row_number, 'chembl_id': cid, 'smiles': smiles, 'errors': errors}, ensure_ascii=False) + '\n')
    for name in files:
        temp[name].replace(files[name])
    stats = {'schema': 1, 'input_sha256': digest, 'records': len(rows), 'safe_accepted': safe_count,
             'npe_accepted': npe_count, 'atom_counts': dict(atoms), 'ring_counts': dict(rings)}
    stat_temp = Path(str(stat_path) + '.tmp')
    stat_temp.write_text(json.dumps(stats), encoding='utf-8')
    stat_temp.replace(stat_path)
    return number, stats

def main():
    CACHE.mkdir(parents=True, exist_ok=True)
    source = OUTPUTS / 'chembl_37_smiles.csv'
    total = json.loads((OUTPUTS / 'dataset_metadata.json').read_text())['exported_rows']
    with source.open('rb') as f:
        source_sha256 = hashlib.file_digest(f, 'sha256').hexdigest()
    counts = Counter()
    atom_counts, ring_counts = Counter(), Counter()
    completed = set()
    started = time.time()
    workers = 12

    def collect(future):
        number, stats = future.result()
        completed.add(number)
        counts.update({key: stats[key] for key in ['records', 'safe_accepted', 'npe_accepted']})
        atom_counts.update(stats['atom_counts'])
        ring_counts.update(stats['ring_counts'])
        progress = {'stage': 'full_corpus_preparation', 'total_records': total, **dict(counts),
                    'completed_shards': len(completed), 'elapsed_seconds': time.time() - started}
        (CACHE / 'progress.json').write_text(json.dumps(progress, indent=2), encoding='utf-8')
        if len(completed) % 10 == 0:
            print(json.dumps(progress), flush=True)

    with ProcessPoolExecutor(max_workers=workers) as pool, source.open(encoding='utf-8', newline='') as f:
        reader = csv.DictReader(f)
        pending = set()
        batch = []
        number = 0
        for index, row in enumerate(reader, 1):
            batch.append((index, row['chembl_id'], row['canonical_smiles']))
            if len(batch) == 5000:
                pending.add(pool.submit(prepare_batch, (number, batch)))
                number += 1
                batch = []
                if len(pending) >= workers * 2:
                    done, pending = wait(pending, return_when=FIRST_COMPLETED)
                    for future in done:
                        collect(future)
        if batch:
            pending.add(pool.submit(prepare_batch, (number, batch)))
            number += 1
        while pending:
            done, pending = wait(pending, return_when=FIRST_COMPLETED)
            for future in done:
                collect(future)
    assert counts['records'] == total
    # Deterministic ordering for tied frequencies, independent of worker completion order.
    sorted_rings = sorted(ring_counts.items(), key=lambda item: (-item[1], item[0]))
    result = {'schema': 1, 'source': str(source), 'source_sha256': source_sha256, 'total_records': total,
              **dict(counts), 'shards': number, 'workers': workers, 'elapsed_seconds': time.time() - started,
              'atom_counts': dict(sorted(atom_counts.items())), 'ring_counts': sorted_rings,
              'safe_rejected': total - counts['safe_accepted'], 'npe_rejected': total - counts['npe_accepted']}
    (CACHE / 'census.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
    progress = dict(result)
    progress.pop('ring_counts')
    progress['stage'] = 'preparation_complete'
    (CACHE / 'progress.json').write_text(json.dumps(progress, indent=2), encoding='utf-8')
    print(json.dumps(progress), flush=True)

if __name__ == '__main__':
    mp.freeze_support()
    main()
