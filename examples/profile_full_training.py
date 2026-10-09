"""Measure full-training feasibility without loading the full NPE graph corpus.

The subsets are resource probes, NOT final scientific training runs. Each case
runs in a fresh process; RSS includes descendants (shared pages may be counted
more than once). The supervisor stops its own process tree before memory is low.
Install the optional analysis dependencies before running this script.
"""
import argparse
from dataclasses import asdict
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import random
import subprocess
import sys
import time

import psutil

GIB = 1024 ** 3


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')


def child(args):
    from rdkit import RDLogger
    from molecular_tokenizer import MolecularTokenizer
    from molecular_tokenizer._vendor import graphbpe
    RDLogger.DisableLog('rdApp.*')
    target = args.output
    events = (target / 'events.jsonl').open('w', encoding='utf-8', buffering=1)
    started = time.perf_counter()

    def event(kind, **fields):
        events.write(json.dumps(dict(event=kind, seconds=time.perf_counter()-started, **fields))+'\n')

    # Instrument progress only; retain the installed training implementation.
    original_tqdm = graphbpe.tqdm

    class TimedProgress(original_tqdm):
        def __init__(self, *a, **kw):
            super().__init__(*a, **kw)
            self.phase = kw.get('desc', '')
            event('phase', phase=self.phase, total=self.total)

        def update(self, n=1):
            result = super().update(n)
            if getattr(self, 'phase', '') == 'Training for subgraph vocabulary':
                event('merge', learned=self.n, total=self.total)
            return result

    graphbpe.tqdm = TimedProgress
    source = args.input.read_text(encoding='utf-8').splitlines()
    event('loaded', rows=len(source), rss=psutil.Process().memory_info().rss)
    rejected = []
    model = MolecularTokenizer(fragmentation=args.fragmentation, representation='safe', bpe_scope='fragment')
    options = dict(num_workers=args.workers, skip_invalid=True,
                   on_reject_with_index=lambda s, r, i: rejected.append(dict(row=i, reason=r)))
    if args.fragmentation == 'npe':
        options.update(motif_vocab_size=args.motifs, ring_vocab_size=args.rings)
        if args.storage == 'sqlite':
            options.update(npe_storage='sqlite',work_dir=target/'checkpoint',
                           progress=lambda record:event('disk_progress',detail=record))
    train_started = time.perf_counter()
    report = model.train(source, args.vocab_size, **options)
    train_seconds = time.perf_counter()-train_started
    model.save(target / 'model')
    event('trained', report=asdict(report), train_seconds=train_seconds)
    # Tiny in-sample smoke test only, never presented as an accuracy benchmark.
    successes = 0
    errors = []
    for smi in source[:100]:
        try:
            model.encode(smi)
            successes += 1
        except Exception as exc:
            errors.append(str(exc))
    write_json(target / 'result.json', dict(report=asdict(report), train_seconds=train_seconds,
               total_seconds=time.perf_counter()-started, rejected=rejected,
               smoke_test=dict(scope='first 100 training rows; not held out', passed=successes, errors=errors)))
    event('finished')
    events.close()


def stop_tree(process):
    try:
        children = psutil.Process(process.pid).children(recursive=True)
    except psutil.NoSuchProcess:
        children = []
    for item in reversed(children):
        try:
            item.kill()
        except psutil.NoSuchProcess:
            pass
    if process.poll() is None:
        process.kill()
    process.wait()
    psutil.wait_procs(children, timeout=5)


def supervise(command, directory, args):
    peak_rss = peak_parent = 0
    reason = None
    started = time.perf_counter()
    with (directory / 'stdout.log').open('w', encoding='utf-8') as stdout, \
         (directory / 'stderr.log').open('w', encoding='utf-8') as stderr, \
         (directory / 'resources.jsonl').open('w', encoding='utf-8', buffering=1) as metrics:
        env = dict(os.environ, PYTHONUNBUFFERED='1', PYTHONIOENCODING='utf-8',
                   TOKENIZERS_PARALLELISM='false', RAYON_NUM_THREADS=str(args.workers))
        process = subprocess.Popen(command, stdout=stdout, stderr=stderr, env=env)
        try:
            while process.poll() is None:
                rss = parent_rss = 0
                try:
                    root = psutil.Process(process.pid)
                    parent_rss = root.memory_info().rss
                    for item in [root]+root.children(recursive=True):
                        try:
                            rss += item.memory_info().rss
                        except psutil.NoSuchProcess:
                            pass
                except psutil.NoSuchProcess:
                    pass
                elapsed = time.perf_counter()-started
                available = psutil.virtual_memory().available
                peak_rss = max(peak_rss, rss)
                peak_parent = max(peak_parent, parent_rss)
                metrics.write(json.dumps(dict(seconds=elapsed, tree_rss=rss,
                    parent_rss=parent_rss, available=available))+'\n')
                if rss > args.max_rss_gib*GIB:
                    reason = 'process-tree RSS limit'
                elif available < args.reserve_gib*GIB:
                    reason = 'available-memory reserve'
                elif elapsed > args.timeout:
                    reason = 'case time limit'
                elif psutil.disk_usage(str(directory)).free < getattr(args, 'min_free_disk_gib', 10)*GIB:
                    reason = 'free-disk reserve'
                if reason:
                    stop_tree(process)
                    break
                time.sleep(0.5)
        finally:
            if process.poll() is None:
                stop_tree(process)
    result = dict(exit_code=process.returncode, stop_reason=reason,
                  peak_tree_rss_gib=peak_rss/GIB, peak_parent_rss_gib=peak_parent/GIB,
                  elapsed_seconds=time.perf_counter()-started)
    write_json(directory/'resources_summary.json', result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--metadata', type=Path)
    parser.add_argument('--workers', type=int, default=4)
    parser.add_argument('--seed', type=int, default=20261009)
    parser.add_argument('--max-rss-gib', type=float, default=10)
    parser.add_argument('--reserve-gib', type=float, default=5)
    parser.add_argument('--timeout', type=float, default=1200)
    parser.add_argument('--cases', default='npe:5000:350,npe:20000:350,npe:50000:350,npe:5000:3000,brics:20000:0')
    parser.add_argument('--child', action='store_true')
    parser.add_argument('--fragmentation', choices=['brics', 'npe'])
    parser.add_argument('--motifs', type=int, default=350)
    parser.add_argument('--rings', type=int, default=300)
    parser.add_argument('--vocab-size', type=int, default=3000)
    parser.add_argument('--storage', choices=['memory','sqlite'], default='memory')
    args = parser.parse_args()
    if args.workers < 1 or args.timeout <= 0 or args.reserve_gib <= 0 or args.max_rss_gib <= 0:
        parser.error('Resource limits and workers must be positive')
    if args.child:
        child(args)
        return
    if args.metadata is None:
        parser.error('--metadata is required for the verified ChEMBL corpus')
    if args.output.exists():
        parser.error('Output directory must not exist')
    source_info = json.loads(args.metadata.read_text(encoding='utf-8'))
    with args.input.open('rb') as stream:
        checksum = hashlib.file_digest(stream, 'sha256').hexdigest()
    if checksum != source_info['outputs']['chembl_37_smiles.smi']['sha256']:
        raise ValueError('Input SHA256 does not match the ChEMBL manifest')
    args.output.mkdir(parents=True)
    count = source_info['exported_rows']
    cases = []
    for item in args.cases.split(','):
        frag, size, motifs = item.split(':')
        size, motifs = int(size), int(motifs)
        if frag not in {'npe', 'brics'} or not 1 <= size <= count:
            raise ValueError(f'Invalid resource-probe case: {item}')
        cases.append(dict(fragmentation=frag, rows=size, motifs=motifs))
    # Nested uniform samples from the full source; restore original row order.
    sample_indices = random.Random(args.seed).sample(range(count), max(c['rows'] for c in cases))
    wanted = set(sample_indices)
    selected = {}
    actual_rows = 0
    with args.input.open(encoding='utf-8') as source:
        for i, line in enumerate(source):
            actual_rows += 1
            if i in wanted:
                selected[i] = line.rstrip('\r\n')
    if actual_rows != count:
        raise ValueError('Row count does not match manifest')
    manifest = dict(scope='Resource assessment only; no final full-corpus model',
                    created_utc=datetime.now(timezone.utc).isoformat(), storage=args.storage, source_rows=count,
                    source_sha256=checksum, seed=args.seed, workers=args.workers,
                    sample_method='nested uniform sample without replacement; restored to source order',
                    hardware=dict(platform=platform.platform(), logical_cpus=psutil.cpu_count(),
                        physical_cpus=psutil.cpu_count(logical=False), total_ram_gib=psutil.virtual_memory().total/GIB,
                        available_ram_gib=psutil.virtual_memory().available/GIB,
                        free_disk_gib=psutil.disk_usage(str(args.output)).free/GIB),
                    limits=dict(max_tree_rss_gib=args.max_rss_gib, reserve_gib=args.reserve_gib,
                                timeout_seconds=args.timeout), cases=[])
    write_json(args.output/'manifest.json', manifest)
    for case in cases:
        name = f"{case['fragmentation']}_{case['rows']}_motifs{case['motifs']}"
        directory = args.output/name
        directory.mkdir()
        indices = sorted(sample_indices[:case['rows']])
        sample_file = directory/'input.smi'
        sample_file.write_text(''.join(selected[i]+'\n' for i in indices), encoding='utf-8')
        (directory/'source_rows.txt').write_text(''.join(str(i+1)+'\n' for i in indices), encoding='utf-8')
        command = [sys.executable, str(Path(__file__).resolve()), '--child', '--input', str(sample_file),
                   '--output', str(directory), '--fragmentation', case['fragmentation'],
                   '--workers', str(args.workers), '--motifs', str(case['motifs']),
                   '--rings', str(args.rings), '--vocab-size', str(args.vocab_size), '--storage',args.storage]
        print(f'Starting resource probe {name}', flush=True)
        resources = supervise(command, directory, args)
        result = dict(case, directory=name, resources=resources)
        if (directory/'result.json').exists():
            result['training'] = json.loads((directory/'result.json').read_text(encoding='utf-8'))
        manifest['cases'].append(result)
        write_json(args.output/'manifest.json', manifest)
        print(json.dumps(result), flush=True)


if __name__ == '__main__':
    main()
