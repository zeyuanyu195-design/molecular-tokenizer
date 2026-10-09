"""Stream a verified ChEMBL corpus through frozen tokenizers with fidelity audits.

Example (all rows, in-corpus description, not a generalization benchmark):
python examples/analyze_corpus.py --input ../chembl_37_smiles.smi \
  --metadata ../dataset_metadata.json --model brics_safe=trained_models/full_brics \
  --model npe_safe=trained_models/full_npe --output analysis_runs/full --workers 4

Without --model, this produces a full-corpus chemistry census. Raw row metrics
are local gzip CSV; aggregate JSON and vocabulary CSV can be published.
"""
import argparse
from collections import Counter
import csv
from dataclasses import asdict
from datetime import datetime, timezone
import gzip
import hashlib
import itertools
import json
import math
import multiprocessing as mp
from pathlib import Path
import time

from rdkit import Chem, RDLogger
from molecular_tokenizer import MolecularTokenizer

MODELS = {}


def initialize(specs):
    RDLogger.DisableLog('rdApp.*')
    global MODELS
    MODELS = {name: MolecularTokenizer.load(path) for name, path in specs}
    for model in MODELS.values():
        # Explicitly compare decoded molecules below, including mismatches.
        # This is not permission to count lossy encodings as successful.
        model.strict = False


def analyze_row(item):
    row, smiles = item
    result = dict(row=row, smiles_characters=len(smiles), models={})
    mol = Chem.MolFromSmiles(smiles) if smiles and smiles == smiles.strip() else None
    if mol is None:
        result['input_valid'] = False
        return result
    canonical = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    flags = [name for name, condition in (
        ('isotope', any(a.GetIsotope() for a in mol.GetAtoms())),
        ('charge', any(a.GetFormalCharge() for a in mol.GetAtoms())),
        ('stereo', any(a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED for a in mol.GetAtoms())
         or any(b.GetStereo() != Chem.BondStereo.STEREONONE for b in mol.GetBonds())),
        ('disconnected', len(Chem.GetMolFrags(mol)) > 1),
    ) if condition]
    result.update(input_valid=True, atoms=mol.GetNumAtoms(), heavy_atoms=mol.GetNumHeavyAtoms(),
                  components=len(Chem.GetMolFrags(mol)), flags=flags)
    for name, model in MODELS.items():
        value = {'status': 'error'}
        try:
            encoded = model.encode(smiles)
            edges = [[e.source, e.target, e.bond_type, e.source_attachment, e.target_attachment]
                     for e in encoded.connections]
            payload = json.dumps({'ids': list(encoded.token_ids), 'edges': edges}, separators=(',', ':'))
            value.update(tokens=len(encoded.token_ids), connections=len(edges),
                         payload_bytes=len(payload.encode('utf-8')),
                         scalar_slots=len(encoded.token_ids)+5*len(edges))
            decoded = model.decode(encoded)
            reconstructed = Chem.MolFromSmiles(decoded)
            exact = reconstructed is not None and Chem.MolToSmiles(
                reconstructed, canonical=True, isomericSmiles=True) == canonical
            value['status'] = 'exact' if exact else 'mismatch'
            if exact:
                value['ids'] = encoded.token_ids
                if model.representation == 'safe':
                    text = model.backend.tokenizer.decode(list(encoded.token_ids), skip_special_tokens=False)
                    value['safe_blocks'] = len(text.split('.'))
            else:
                value['reason'] = 'canonical_isomeric_smiles_mismatch'
        except Exception as error:
            value['reason'] = type(error).__name__
            value['detail'] = str(error)[:300]
        result['models'][name] = value
    return result


class Aggregate:
    def __init__(self, names):
        self.rows = self.valid = 0
        self.atom_hist = Counter()
        self.char_hist = Counter()
        self.groups = Counter()
        self.models = {name: dict(status=Counter(), reasons=Counter(), tokens=Counter(),
            connections=Counter(), payload=Counter(), ratio_tenths=Counter(), usage=Counter(),
            joint=Counter(), groups={}, total_atoms=0, total_tokens=0, sum_atom_token_ratio=0.,
            examples=[]) for name in names}
        self.paired = {name: dict(n=0, tokens=0, connections=0, payload_bytes=0) for name in names}

    def add(self, row):
        self.rows += 1
        self.char_hist[row['smiles_characters']] += 1
        if not row['input_valid']:
            for value in self.models.values():
                value['status']['invalid_input'] += 1
            return
        self.valid += 1
        self.atom_hist[row['atoms']] += 1
        self.groups.update(row['flags'])
        all_exact = bool(self.models) and all(v['status'] == 'exact' for v in row['models'].values())
        for name, value in row['models'].items():
            stats = self.models[name]
            stats['status'][value['status']] += 1
            for group in row['flags']:
                stats['groups'].setdefault(group, Counter())[value['status']] += 1
            if value['status'] != 'exact':
                stats['reasons'][value.get('reason', 'unknown')] += 1
                if len(stats['examples']) < 12:
                    stats['examples'].append(dict(source_row=row['row'], **value))
                continue
            stats['tokens'][value['tokens']] += 1
            stats['connections'][value['connections']] += 1
            stats['payload'][value['payload_bytes']] += 1
            stats['usage'].update(value['ids'])
            ratio = row['atoms']/value['tokens']
            stats['ratio_tenths'][int(ratio*10)] += 1
            # Binned, not sampled scatter: every successful molecule contributes.
            stats['joint'][(row['atoms']//5, value['tokens']//5)] += 1
            stats['total_atoms'] += row['atoms']
            stats['total_tokens'] += value['tokens']
            stats['sum_atom_token_ratio'] += ratio
            if all_exact:
                pair = self.paired[name]
                pair['n'] += 1
                for key in ('tokens', 'connections', 'payload_bytes'):
                    pair[key] += value[key]

    def result(self):
        output = dict(rows=self.rows, valid_molecules=self.valid, invalid_molecules=self.rows-self.valid,
                      atom_histogram=dict(self.atom_hist), smiles_length_histogram=dict(self.char_hist),
                      chemistry_groups=dict(self.groups), models={})
        for name, value in self.models.items():
            n = value['status']['exact']
            total = value['total_tokens']
            counts = sorted(value['usage'].values(), reverse=True)
            entropy = -sum(c/total*math.log2(c/total) for c in counts) if total else None
            output['models'][name] = dict(
                status=dict(value['status']), reasons=dict(value['reasons']),
                exact_fraction_all_rows=n/self.rows if self.rows else None,
                exact_fraction_valid_inputs=n/self.valid if self.valid else None,
                token_histogram=dict(value['tokens']), connection_histogram=dict(value['connections']),
                payload_bytes_histogram=dict(value['payload']),
                atom_token_ratio_histogram_tenths=dict(value['ratio_tenths']),
                atom_token_joint_bins5=[[a, t, c] for (a,t),c in value['joint'].items()],
                token_usage=dict(value['usage']), usage_entropy_bits=entropy,
                groups={k:dict(v) for k,v in value['groups'].items()},
                mean_tokens_on_exact=total/n if n else None,
                mean_atom_token_ratio_on_exact=value['sum_atom_token_ratio']/n if n else None,
                ratio_of_total_atoms_to_total_tokens=value['total_atoms']/total if total else None,
                common_exact_subset=self.paired[name], failure_examples=value['examples'])
        return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--metadata', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--model', action='append', default=[], help='LABEL=MODEL_DIRECTORY')
    parser.add_argument('--workers', type=int, default=4)
    parser.add_argument('--max-molecules', type=int)
    parser.add_argument('--scope', choices=['training_corpus', 'held_out', 'unknown'], default='training_corpus')
    args = parser.parse_args()
    if args.workers < 1 or (args.max_molecules is not None and args.max_molecules < 1):
        parser.error('Workers and max-molecules must be positive')
    if args.output.exists():
        parser.error('Output directory must not exist')
    specs = [item.split('=', 1) for item in args.model]
    if any(len(s) != 2 for s in specs) or len({s[0] for s in specs}) != len(specs):
        parser.error('Model labels must be unique LABEL=PATH pairs')
    source_info = json.loads(args.metadata.read_text(encoding='utf-8'))
    with args.input.open('rb') as stream:
        source_sha = hashlib.file_digest(stream, 'sha256').hexdigest()
    if source_sha != source_info['outputs']['chembl_37_smiles.smi']['sha256']:
        raise ValueError('Input SHA256 differs from source manifest')
    args.output.mkdir(parents=True)
    aggregate = Aggregate([name for name, _ in specs])
    started = time.perf_counter()
    with args.input.open(encoding='utf-8') as source, \
         gzip.open(args.output/'molecules.csv.gz', 'wt', encoding='utf-8', newline='') as destination:
        fields = ['source_row','valid','atoms','smiles_characters','flags','model','status',
                  'tokens','connections','payload_bytes','scalar_slots','safe_blocks','reason']
        writer = csv.DictWriter(destination, fieldnames=fields)
        writer.writeheader()
        rows = ((i, s.rstrip('\r\n')) for i,s in enumerate(itertools.islice(source, args.max_molecules),1))
        with mp.Pool(args.workers, initializer=initialize, initargs=(specs,)) as pool:
            for row in pool.imap(analyze_row, rows, chunksize=32):
                aggregate.add(row)
                base = dict(source_row=row['row'], valid=row['input_valid'], atoms=row.get('atoms'),
                            smiles_characters=row['smiles_characters'], flags='|'.join(row.get('flags', [])))
                if not specs:
                    writer.writerow(base)
                for name, _ in specs:
                    values = row['models'].get(name, {'status':'invalid_input'})
                    writer.writerow(dict(base, model=name, **{k:v for k,v in values.items() if k in fields}))
                if aggregate.rows % 10000 == 0:
                    print(json.dumps(dict(rows=aggregate.rows, seconds=time.perf_counter()-started)), flush=True)
    expected = min(args.max_molecules or source_info['exported_rows'], source_info['exported_rows'])
    if aggregate.rows != expected:
        raise ValueError(f'Expected {expected} rows, read {aggregate.rows}')
    result = aggregate.result()
    result.update(source_sha256=source_sha, source_rows=source_info['exported_rows'],
        selection='all rows' if args.max_molecules is None else f'first {args.max_molecules} source rows',
        scope=args.scope, created_utc=datetime.now(timezone.utc).isoformat(),
        elapsed_seconds=time.perf_counter()-started, workers=args.workers,
        definitions=dict(exact='canonical isomeric SMILES identity after decoding; failed inputs stay in denominator',
            atoms='RDKit GetNumAtoms, including explicit atoms; implicit hydrogens excluded',
            ratio='atoms / SAFE sequence tokens or atoms / DemoDiff nodes; distinct metrics, not byte compression',
            payload='UTF-8 compact JSON {ids:[...],edges:[[source,target,bond_type,source_attachment,target_attachment],...]}; excludes model storage and metadata',
            common_subset='Only molecules reconstructed exactly by every supplied model',
            chemistry_groups='Overlapping RDKit-based subsets; counts must not be summed'))
    for name, path in specs:
        model = MolecularTokenizer.load(path)
        result['models'][name]['configuration'] = dict(fragmentation=model.fragmentation,
            representation=model.representation, vocab_size=model.vocab_size, tokenizer_id=model.tokenizer_id)
        partition = getattr(model.backend,'npe',None)
        result['models'][name]['configuration']['partition_vocab_size'] = partition.vocab_size if partition else None
        result['models'][name]['configuration']['saved_model_bytes'] = sum(
            p.stat().st_size for p in Path(path).iterdir() if p.is_file() and p.name!='training_report.json')
        provenance=Path(path)/'training_report.json'
        if provenance.exists():
            result['models'][name]['training_provenance']=json.loads(provenance.read_text(encoding='utf-8'))
        if model.representation == 'safe':
            units = sorted(model.backend.tokenizer.get_vocab().items(), key=lambda p:p[1])
        else:
            units = [(token,i) for i,token in enumerate(model.backend.engine.vocab_node)]
        # Numeric filenames avoid interpreting caller-supplied labels as paths.
        filename = f'vocabulary_{list(aggregate.models).index(name)}.csv'
        result['models'][name]['vocabulary_file'] = filename
        with (args.output/filename).open('w', encoding='utf-8', newline='') as sink:
            writer = csv.writer(sink)
            writer.writerow(['id','unit','exact_reconstruction_usage','string_characters','motif_atoms'])
            for token, i in units:
                motif_atoms = (model.backend.engine.vocab_node_stats[token][0]
                               if model.representation == 'demodiff' else '')
                writer.writerow([i, token, aggregate.models[name]['usage'][i], len(token), motif_atoms])
    (args.output/'summary.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k in ['rows','valid_molecules','invalid_molecules','elapsed_seconds']}))


if __name__ == '__main__':
    main()
