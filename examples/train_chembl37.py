"""Train separate SAFE/NPE vocabularies on the prepared full ChEMBL 37 corpus."""
import argparse
from collections import Counter
import csv
import hashlib
import json
import multiprocessing as mp
from pathlib import Path
import random
import sys
from rdkit import RDLogger
from time import monotonic

from molecular_tokenizer import MolecularTokenizer

PROJECT = Path(__file__).resolve().parents[1]
OUTPUTS = PROJECT.parent
SOURCE = OUTPUTS / 'chembl_37_smiles.smi'
SOURCE_METADATA = OUTPUTS / 'dataset_metadata.json'
TARGETS = {
    'safe': (PROJECT / 'trained_models' / 'chembl37_safe_vocab3000', 3000,
             {'num_workers': 4}),
    'npe': (PROJECT / 'trained_models' / 'chembl37_npe_vocab3000ring300', 3000,
            {'ring_vocab_size': 300, 'num_workers': 2}),
}

def sha256(path):
    with path.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()

def load_ids():
    with (OUTPUTS / 'chembl_37_smiles.csv').open(encoding='utf-8', newline='') as source_csv:
        return [record['chembl_id'] for record in csv.DictReader(source_csv)]

def read_smiles(progress=False):
    metadata = json.loads(SOURCE_METADATA.read_text(encoding='utf-8'))
    expected_count = metadata['exported_rows']
    expected_sha = metadata['outputs']['chembl_37_smiles.smi']['sha256']
    if SOURCE.stat().st_size != metadata['outputs']['chembl_37_smiles.smi']['bytes'] or sha256(SOURCE) != expected_sha:
        raise RuntimeError('ChEMBL SMI source differs from the verified step 1 dataset')
    started = monotonic()
    count = 0
    source_csv = (OUTPUTS / 'chembl_37_smiles.csv').open(encoding='utf-8', newline='')
    with SOURCE.open(encoding='utf-8') as f:
        records = csv.DictReader(source_csv)
        try:
          for line, record in zip(f, records, strict=True):
            smiles = line.rstrip('\r\n')
            if record['canonical_smiles'] != smiles:
                raise RuntimeError(f'SMILES/ID source mismatch at row {count + 1}')
            if not smiles:
                raise RuntimeError(f'Blank SMILES at source row {count + 1}')
            count += 1
            if progress and count % 100000 == 0:
                print(f'Read {count:,}/{expected_count:,} molecules ({count / expected_count:.1%}); elapsed {(monotonic()-started)/60:.1f} min', flush=True)
            yield {'chembl_id': record['chembl_id'], 'row': count, 'smiles': smiles}
        finally:
            source_csv.close()
    if count != expected_count:
        raise RuntimeError(f'Expected {expected_count:,} SMILES, read {count:,}')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('backend', choices=['safe', 'npe'])
    parser.add_argument('--vocab-size', type=int, default=3000)
    parser.add_argument('--num-workers', type=int)
    parser.add_argument('--max-rows', type=int,
                        help='Use a deterministic reservoir sample of this many rows (for memory-limited NPE runs)')
    parser.add_argument('--seed', type=int, default=42)
    args = parser.parse_args()
    if args.max_rows is not None and args.max_rows < 1:
        parser.error('--max-rows must be positive')
    target, _, default_options = TARGETS[args.backend]
    options = {'num_workers': args.num_workers or default_options['num_workers']}
    if args.backend == 'npe':
        options['ring_vocab_size'] = default_options.get('ring_vocab_size', 300)
    target = target.parent / target.name.replace('3000', str(args.vocab_size))
    if target.exists():
        raise FileExistsError(f'Refusing to overwrite existing training run: {target}')

    verified = json.loads(SOURCE_METADATA.read_text(encoding='utf-8'))
    report = {'backend': args.backend, 'status': 'training', 'source': str(SOURCE),
              'source_rows': verified['exported_rows'], 'source_sha256': verified['outputs']['chembl_37_smiles.smi']['sha256'],
              'vocab_size_target': args.vocab_size, 'backend_options': options}
    print(json.dumps(report, ensure_ascii=False), flush=True)
    rejected_file = target.parent / f'{target.name}_rejected.tsv'
    rejected_stats = Counter()
    rejected_examples = []
    rejected_count = 0
    chembl_ids = load_ids()
    if len(chembl_ids) != verified['exported_rows']:
        raise RuntimeError('ChEMBL ID row count differs from the verified dataset manifest')
    target.parent.mkdir(parents=True, exist_ok=True)
    rejected_sink = rejected_file.open('w', encoding='utf-8', newline='\n')
    rejected_sink.write('source_row\tchembl_id\treason\n')
    training_source_rows = []
    def on_reject(smiles, reason, training_index):
        nonlocal rejected_count
        source_row = training_source_rows[training_index - 1]
        rejected_count += 1
        reason_type = reason.split(':', 1)[0]
        rejected_stats[reason_type] += 1
        if len(rejected_examples) < 25:
            rejected_examples.append({'source_row': source_row, 'chembl_id': chembl_ids[source_row - 1],
                                      'smiles': smiles, 'reason': reason[:1000]})
        chembl_id = chembl_ids[source_row - 1]
        rejected_sink.write(f'{source_row}\t{chembl_id}\t{reason.replace(chr(9), " ").replace(chr(10), " ")[:500]}\n')
        if rejected_count <= 10 or rejected_count % 1000 == 0:
            print(f'Rejected {rejected_count:,} molecules; current ChEMBL ID={chembl_id}, reason={reason_type}', flush=True)
    def selected_records():
        records = read_smiles(progress=True)
        if args.max_rows is None or args.max_rows >= verified['exported_rows']:
            yield from records
            return
        rng = random.Random(args.seed)
        reservoir = []
        for seen, record in enumerate(records):
            if len(reservoir) < args.max_rows:
                reservoir.append(record)
            else:
                replacement = rng.randrange(seen + 1)
                if replacement < args.max_rows:
                    reservoir[replacement] = record
        yield from sorted(reservoir, key=lambda item: item['row'])

    def records_with_id():
        for record in selected_records():
            training_source_rows.append(record['row'])
            yield record['smiles']
    model = MolecularTokenizer(args.backend)
    started = monotonic()
    RDLogger.DisableLog('rdApp.*')
    try:
        training = model.train(records_with_id(), args.vocab_size, skip_invalid=True,
                               on_reject_with_index=on_reject, **options)
    finally:
        rejected_sink.close()
    if not rejected_count:
        rejected_file.unlink()
    model.save(target)

    if args.backend == 'safe':
        tokens = model.backend.tokenizer.get_vocab()
        with (target / 'vocabulary.tsv').open('w', encoding='utf-8', newline='\n') as f:
            f.write('token_id\ttoken\n')
            for token, token_id in sorted(tokens.items(), key=lambda item: item[1]):
                f.write(f'{token_id}\t{token}\n')
    else:
        engine = model.backend.engine
        with (target / 'vocabulary.tsv').open('w', encoding='utf-8', newline='\n') as f:
            f.write('token_id\tmotif\tatom_count\tfrequency\n')
            for token_id, motif in enumerate(engine.vocab_node):
                atom_count, frequency = engine.vocab_node_stats[motif]
                f.write(f'{token_id}\t{motif}\t{atom_count}\t{frequency}\n')
    selected_rows_file = None
    if args.max_rows is not None and args.max_rows < verified['exported_rows']:
        selected_rows_file = target / 'training_selection.tsv'
        with selected_rows_file.open('w', encoding='utf-8', newline='\n') as f:
            f.write('source_row\tchembl_id\n')
            for source_row in training_source_rows:
                f.write(f'{source_row}\t{chembl_ids[source_row - 1]}\n')

    record = {'status': 'complete', 'backend': args.backend, 'source': str(SOURCE),
              'source_rows': verified['exported_rows'], 'source_sha256': verified['outputs']['chembl_37_smiles.smi']['sha256'],
              'training_input_rows': len(training_source_rows),
              'accepted_molecules': training.molecules,
              'vocab_size_target': args.vocab_size, 'actual_node_or_sequence_vocab_size': training.vocab_size,
              'ring_vocab_size': training.ring_vocab_size, 'edge_vocab_size': training.edge_vocab_size,
              'training_seconds': monotonic() - started, 'strict_fidelity': model.strict,
              'tokenizer_id': model.tokenizer_id, 'model_files': sorted(p.name for p in target.iterdir()),
              'training_options': options, 'source_order': 'official ChEMBL 37 export order',
              'randomized_subset': selected_rows_file is not None,
              'sampling_method': 'uniform reservoir sample, restored to source order' if selected_rows_file else 'full source order',
              'sampling_seed': args.seed if selected_rows_file else None,
              'training_selection_file': str(selected_rows_file) if selected_rows_file else None,
              'training_validation_split': False,
              'rejected_molecules': rejected_count, 'rejection_reason_counts': dict(rejected_stats),
              'rejection_examples': rejected_examples,
              'rejection_log': str(rejected_file) if rejected_count else None,
              'warnings': ['Invalid/unconvertible molecules were explicitly logged and skipped.'
                           ] if rejected_count else []}
    (target / 'training_manifest.json').write_text(json.dumps(record, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(record, ensure_ascii=False, indent=2), flush=True)

if __name__ == '__main__':
    mp.freeze_support()
    main()
