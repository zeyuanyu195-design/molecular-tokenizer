"""Exercise the shared API and audit a bounded sample of the prepared ChEMBL dataset."""
import argparse
from dataclasses import asdict
from itertools import islice
import json
import multiprocessing as mp
from pathlib import Path
import tempfile
from time import perf_counter

from molecular_tokenizer import MolecularTokenizer, MoleculeEncoding

PROJECT = Path(__file__).resolve().parents[1]
OUTPUTS = PROJECT.parent

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--limit', type=int, default=100)
    args = parser.parse_args()
    if args.limit < 1:
        parser.error('--limit must be positive')
    with (OUTPUTS / 'chembl_37_smiles.smi').open(encoding='utf-8') as source:
        data = [line.strip() for line in islice(source, args.limit)]
    validation = PROJECT / 'validation'
    validation.mkdir(exist_ok=True)
    run = Path(tempfile.mkdtemp(prefix='demo-', dir=validation))

    safe = MolecularTokenizer('safe')
    safe_report = safe.train(iter(data), vocab_size=256)
    # Import the existing GraphBPE vocabulary into the same facade.
    npe = MolecularTokenizer.from_graphbpe(OUTPUTS / 'DemoDiff/data/tokenizer/vocab3000ring300/pretrain-token')
    npe_training = MolecularTokenizer('npe')
    toy = ['CCO', 'CCCO', 'CCN', 'CCCN', 'CC(=O)O', 'CCOC(=O)C',
           'c1ccccc1', 'Cc1ccccc1', 'Oc1ccccc1', 'c1ccncc1', 'CCCl', 'CCF']
    npe_report = npe_training.train(toy, vocab_size=16, ring_vocab_size=2, num_workers=2)
    npe_training.save(run / 'npe_trained_toy')

    audit = {'sample_size': len(data), 'sample_selection': 'First N records in official ChEMBL 37 export; not a random or full-dataset benchmark',
             'safe_training': asdict(safe_report), 'npe_toy_training': asdict(npe_report), 'backends': {}}
    for name, model in [('safe', safe), ('npe', npe)]:
        path = model.save(run / name)
        loaded = MolecularTokenizer.load(path)
        rows = []
        elapsed = perf_counter()
        for sample in data:
            try:
                payload = model.encode(sample)
                serialized = json.loads(json.dumps(payload.to_dict()))
                restored = MoleculeEncoding.from_dict(serialized)
                decoded = loaded.decode(restored)
                rows.append({'input': sample, 'ok': True, 'token_count': len(payload.token_ids),
                             'connection_count': len(payload.connections), 'decoded': decoded,
                             'encoding': serialized})
            except Exception as error:
                rows.append({'input': sample, 'ok': False, 'error': f'{type(error).__name__}: {error}'})
        successes = sum(row['ok'] for row in rows)
        audit['backends'][name] = {'accepted': successes, 'rejected': len(rows) - successes,
                                  'seconds': perf_counter() - elapsed,
                                  'vocab_size': model.vocab_size, 'model_path': str(path)}
        (run / f'{name}_audit.json').write_text(json.dumps(rows, ensure_ascii=False, indent=2), encoding='utf-8')
    (validation / 'latest_demo_summary.json').write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps(audit, ensure_ascii=False, indent=2))

if __name__ == '__main__':
    mp.freeze_support()
    main()
