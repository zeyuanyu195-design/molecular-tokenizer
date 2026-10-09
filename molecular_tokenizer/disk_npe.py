"""Exact, disk-backed implementation of the installed NPE training loop.

SQLite stores graph states, candidate counts and the inverted index. Only a
bounded batch of graphs is in Python memory. Frequency ties follow the legacy
dictionary insertion order, including deletion and later reactivation. A merge
removes ALL old counts before adding ANY replacement counts, just as in memory.

Checkpoints contain Python/RDKit pickles and must be trusted local files created
by this implementation. SQLite transactions make batch/cursor updates atomic.
"""
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
import multiprocessing as mp
import os
from pathlib import Path
import pickle
import sqlite3
import time
import zlib

from rdkit import Chem, RDLogger, rdBase

from .backends import canonical_molecule
from .types import TokenizerError, UnsupportedMoleculeError
from ._vendor.graphbpe import read_single_molecule, count_subgraph_frequency
from ._vendor.graphbpe_utils import ELEMENTS, count_atom, get_indexed_smiles


def _init_worker():
    RDLogger.DisableLog('rdApp.*')


def _pack(value):
    return zlib.compress(pickle.dumps(value, protocol=5), level=1)


def _unpack(value):
    return pickle.loads(zlib.decompress(value))


def _read_graph(item):
    row, smiles = item
    try:
        canonical_molecule(smiles)
        if '*' in smiles:
            raise UnsupportedMoleculeError('Open dummy attachment points are unsupported')
        mol = Chem.MolFromSmiles(smiles)
        Chem.Kekulize(mol, clearAromaticFlags=True)
        processor, rings = read_single_molecule(smiles, True, False)
        if processor is None:
            raise UnsupportedMoleculeError('NPE graph initialization failed')
        # Legacy train_node transports the initialized processor through a
        # multiprocessing queue before candidate counting. RDKit pickling drops
        # computed properties that affect fragment SMILES for some chiral atoms.
        # Reproduce that boundary rather than silently changing merge counts.
        processor = pickle.loads(pickle.dumps(processor, protocol=5))
        counts, processor = count_subgraph_frequency(processor)
        atoms = Counter(a.GetSymbol() for a in processor.molecule.GetAtoms())
        return row, _pack(processor), counts, dict(atoms), rings, None
    except Exception as error:
        return row, None, None, None, None, f'{type(error).__name__}: {error}'


def _merge_graph(item):
    row, blob, selected = item
    processor = _unpack(blob)
    processor.merge(selected)
    processor = pickle.loads(pickle.dumps(processor, protocol=5))
    frequencies, processor = count_subgraph_frequency(processor)
    return row, _pack(processor), frequencies


class DiskNPETrainer:
    """One local checkpoint; workers never write to the database."""
    def __init__(self, directory, vocab_size, ring_size, *, workers=4, batch_size=256,
                 resume=False, skip_invalid=False, progress=None):
        if workers < 1 or batch_size < 1 or vocab_size < 1 or not 0 <= ring_size <= vocab_size:
            raise TokenizerError('Invalid disk NPE configuration')
        self.directory = Path(directory)
        self.directory.mkdir(parents=True, exist_ok=True)
        filename = self.directory/'npe.sqlite3'
        if filename.exists() and not resume:
            raise TokenizerError('NPE checkpoint already exists; explicitly set resume=True')
        self.lock = (self.directory/'writer.lock').open('a+b')
        self.lock.seek(0,2)
        if self.lock.tell()==0:
            self.lock.write(b'0')
            self.lock.flush()
        self.lock.seek(0)
        try:
            if os.name=='nt':
                import msvcrt
                msvcrt.locking(self.lock.fileno(),msvcrt.LK_NBLCK,1)
            else:
                import fcntl
                fcntl.flock(self.lock.fileno(),fcntl.LOCK_EX|fcntl.LOCK_NB)
        except OSError:
            self.lock.close()
            raise TokenizerError('Another process is using this NPE checkpoint')
        self.workers, self.batch_size = workers, batch_size
        self.skip_invalid, self.progress = skip_invalid, progress
        self.db = sqlite3.connect(filename, timeout=120)
        self.db.execute('PRAGMA journal_mode=WAL')
        self.db.execute('PRAGMA synchronous=FULL')
        self.db.execute('PRAGMA cache_size=-131072')
        self.db.execute('PRAGMA temp_store=FILE')
        self.db.execute('PRAGMA mmap_size=0')
        self.db.executescript('''
            CREATE TABLE IF NOT EXISTS meta(key TEXT PRIMARY KEY, value TEXT NOT NULL);
            CREATE TABLE IF NOT EXISTS inputs(id INTEGER PRIMARY KEY, smiles TEXT NOT NULL);
            CREATE TABLE IF NOT EXISTS rejected(id INTEGER PRIMARY KEY, reason TEXT NOT NULL);
            CREATE TABLE IF NOT EXISTS molecules(id INTEGER PRIMARY KEY, accepted_index INTEGER UNIQUE,
                state BLOB NOT NULL, frequencies BLOB NOT NULL);
            CREATE TABLE IF NOT EXISTS candidates(id INTEGER PRIMARY KEY, smiles TEXT UNIQUE NOT NULL,
                total INTEGER NOT NULL CHECK(total>=0), rank INTEGER NOT NULL);
            CREATE INDEX IF NOT EXISTS choose_candidate ON candidates(total DESC, rank) WHERE total>0;
            CREATE TABLE IF NOT EXISTS occurrences(candidate INTEGER, molecule INTEGER, count INTEGER NOT NULL,
                PRIMARY KEY(candidate,molecule)) WITHOUT ROWID;
            CREATE TABLE IF NOT EXISTS affected(id INTEGER PRIMARY KEY);
            CREATE TABLE IF NOT EXISTS rings(id INTEGER PRIMARY KEY, smiles TEXT UNIQUE, count INTEGER NOT NULL);
            CREATE TABLE IF NOT EXISTS atoms(symbol TEXT PRIMARY KEY, count INTEGER NOT NULL);
            CREATE TABLE IF NOT EXISTS vocab(position INTEGER PRIMARY KEY, smiles TEXT UNIQUE,
                atoms INTEGER NOT NULL, frequency INTEGER NOT NULL);
        ''')
        config = dict(algorithm='disk-npe-exact-v1', vocab_size=vocab_size, ring_size=ring_size,
                      rdkit=rdBase.rdkitVersion, skip_invalid=skip_invalid)
        previous = self.get('config')
        if previous is not None and previous != config:
            self.close()
            raise TokenizerError('Checkpoint configuration/RDKit differs from requested training')
        with self.db:
            self.set('config', config)
            if self.get('stage') is None:
                self.set('stage', 'input')
                self.set('rank', 0)
                self.set('cursor', 0)
                self.set('accepted', 0)
                self.set('merges', 0)

    def close(self):
        if self.db is not None:
            self.db.close()
            self.db = None
        if getattr(self,'lock',None) is not None:
            self.lock.close()
            self.lock = None

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

    def get(self, key, default=None):
        record = self.db.execute('SELECT value FROM meta WHERE key=?',(key,)).fetchone()
        return json.loads(record[0]) if record else default

    def set(self, key, value):
        self.db.execute('INSERT INTO meta VALUES(?,?) ON CONFLICT(key) DO UPDATE SET value=excluded.value',
                        (key,json.dumps(value)))

    def notify(self, event, **values):
        record = dict(event=event, utc=datetime.now(timezone.utc).isoformat(), **values)
        # No source molecule strings in operational progress logs.
        with (self.directory/'progress.jsonl').open('a',encoding='utf-8') as stream:
            stream.write(json.dumps(record)+'\n')
        if self.progress:
            self.progress(record)

    def spool(self, smiles):
        """Hash every source row on resume BEFORE changing graph-training state."""
        digest, count = hashlib.sha256(), 0
        fresh = self.get('stage') == 'input'
        try:
            if fresh:
                self.db.execute('DELETE FROM inputs')
            batch=[]
            for count, value in enumerate(smiles,1):
                if not isinstance(value,str):
                    raise TokenizerError('Disk training requires an iterable of SMILES strings')
                data=value.encode('utf-8')
                digest.update(len(data).to_bytes(8,'big'))
                digest.update(data)
                if fresh:
                    batch.append((count,value))
                    if len(batch)>=4096:
                        self.db.executemany('INSERT INTO inputs VALUES(?,?)',batch)
                        batch.clear()
            identity=dict(rows=count,sha256_length_prefixed_utf8=digest.hexdigest())
            if not count:
                raise TokenizerError('Training corpus is empty')
            if fresh:
                self.db.executemany('INSERT INTO inputs VALUES(?,?)',batch)
                self.set('input',identity)
                self.set('stage','initialize')
                self.db.commit()
            elif identity != self.get('input'):
                raise TokenizerError('Checkpoint source rows/order/contents differ from this training input')
        except BaseException:
            self.db.rollback()
            raise

    def _add_counts(self, graph_frequencies):
        """Positive increments only; first-activation order matches Python dict."""
        totals=Counter()
        for _, frequencies in graph_frequencies:
            totals.update(frequencies)
        rank=self.get('rank')
        identifiers={}
        for smiles,count in totals.items():
            record=self.db.execute('SELECT id,total FROM candidates WHERE smiles=?',(smiles,)).fetchone()
            if record is None:
                rank+=1
                result=self.db.execute('INSERT INTO candidates(smiles,total,rank) VALUES(?,?,?)',(smiles,count,rank))
                identifiers[smiles]=result.lastrowid
            else:
                candidate,previous=record
                if previous==0:
                    rank+=1
                    self.db.execute('UPDATE candidates SET total=?,rank=? WHERE id=?',(count,rank,candidate))
                else:
                    self.db.execute('UPDATE candidates SET total=total+? WHERE id=?',(count,candidate))
                identifiers[smiles]=candidate
        self.set('rank',rank)
        mapped={}
        occurrences=[]
        for row,frequencies in graph_frequencies:
            pairs=[(identifiers[s],count) for s,count in frequencies.items()]
            mapped[row]=_pack(pairs)
            occurrences.extend((candidate,row,count) for candidate,count in pairs)
        self.db.executemany('INSERT INTO occurrences VALUES(?,?,?)',occurrences)
        return mapped

    def initialize_graphs(self,pool):
        cursor=self.get('cursor')
        while True:
            batch=self.db.execute('SELECT id,smiles FROM inputs WHERE id>? ORDER BY id LIMIT ?',
                                  (cursor,self.batch_size)).fetchall()
            if not batch:
                break
            results=pool.map(_read_graph,batch)
            good=[r for r in results if r[-1] is None]
            with self.db:
                accepted=self.get('accepted')
                mapped=self._add_counts([(r[0],r[2]) for r in good])
                for row,blob,freq,atoms,rings,error in results:
                    if error:
                        if not self.skip_invalid:
                            raise TokenizerError(error)
                        self.db.execute('INSERT INTO rejected VALUES(?,?)',(row,error))
                        continue
                    accepted+=1
                    self.db.execute('INSERT INTO molecules VALUES(?,?,?,?)',(row,accepted,blob,mapped[row]))
                    self.db.executemany('INSERT INTO atoms VALUES(?,?) ON CONFLICT(symbol) DO UPDATE SET count=count+excluded.count',atoms.items())
                    self.db.executemany('INSERT INTO rings(smiles,count) VALUES(?,?) ON CONFLICT(smiles) DO UPDATE SET count=count+excluded.count',rings.items())
                cursor=batch[-1][0]
                self.set('accepted',accepted)
                self.set('cursor',cursor)
            if cursor% (self.batch_size*40)==0:
                self.notify('initialize',source_rows=cursor,accepted=accepted,total=self.get('input')['rows'])
        if not self.get('accepted'):
            raise TokenizerError('Training corpus has no accepted molecules')
        atom_counts=dict(self.db.execute('SELECT symbol,count FROM atoms'))
        rings=self.db.execute('SELECT smiles,count FROM rings ORDER BY count DESC,id LIMIT ?',
                              (self.get('config')['ring_size'],)).fetchall()
        stats={symbol:[1,atom_counts.get(symbol,0)] for symbol in ELEMENTS}
        for smiles,count in rings:
            stats[smiles]=[count_atom(smiles),count]
        vocab=[s for s in ELEMENTS+[r[0] for r in rings] if stats[s][1]>0]
        if len(vocab)>self.get('config')['vocab_size']:
            raise TokenizerError('Motif budget is smaller than the initial atom/ring vocabulary')
        with self.db:
            self.db.executemany('INSERT INTO vocab VALUES(?,?,?,?)',
                               [(i,s,*stats[s]) for i,s in enumerate(vocab)])
            self.set('initial_rings',[r[0] for r in rings])
            self.set('stage','choose')
            self.set('cursor',0)
        self.notify('initialized',accepted=self.get('accepted'),initial_vocab=len(vocab))

    def choose(self):
        size=self.db.execute('SELECT count(*) FROM vocab').fetchone()[0]
        selected=self.db.execute('SELECT id,smiles,total FROM candidates WHERE total>0 ORDER BY total DESC,rank LIMIT 1').fetchone()
        with self.db:
            if size>=self.get('config')['vocab_size'] or selected is None:
                self.set('stage','complete')
                return
            candidate,smiles,frequency=selected
            self.db.execute('DELETE FROM affected')
            self.db.execute('INSERT INTO affected SELECT molecule FROM occurrences WHERE candidate=?',(candidate,))
            self.set('selected',dict(id=candidate,smiles=smiles,frequency=frequency))
            self.set('stage','remove')
            self.set('cursor',0)

    def remove_counts(self):
        cursor=self.get('cursor')
        batches=0
        while True:
            batch=self.db.execute('SELECT m.id,m.frequencies FROM affected a JOIN molecules m ON m.id=a.id WHERE a.id>? ORDER BY a.id LIMIT ?',
                                  (cursor,self.batch_size)).fetchall()
            if not batch:
                break
            totals=Counter()
            deletions=[]
            for row,blob in batch:
                for candidate,count in _unpack(blob):
                    totals[candidate]+=count
                    deletions.append((candidate,row))
            with self.db:
                self.db.executemany('UPDATE candidates SET total=total-? WHERE id=?',[(v,k) for k,v in totals.items()])
                self.db.executemany('DELETE FROM occurrences WHERE candidate=? AND molecule=?',deletions)
                cursor=batch[-1][0]
                self.set('cursor',cursor)
            batches+=1
            if batches%40==0:
                self.notify('checkpoint',phase='remove',last_source_row=cursor,merge=self.get('merges')+1)
        with self.db:
            self.set('stage','update')
            self.set('cursor',0)

    def update_graphs(self,pool):
        cursor=self.get('cursor')
        selected=self.get('selected')
        batches=0
        while True:
            batch=self.db.execute('SELECT m.id,m.state FROM affected a JOIN molecules m ON m.id=a.id WHERE a.id>? ORDER BY a.id LIMIT ?',
                                  (cursor,self.batch_size)).fetchall()
            if not batch:
                break
            results=pool.map(_merge_graph,[(row,blob,selected['smiles']) for row,blob in batch])
            with self.db:
                mapped=self._add_counts([(row,freq) for row,_,freq in results])
                self.db.executemany('UPDATE molecules SET state=?,frequencies=? WHERE id=?',
                                    [(blob,mapped[row],row) for row,blob,_ in results])
                cursor=batch[-1][0]
                self.set('cursor',cursor)
            batches+=1
            if batches%40==0:
                self.notify('checkpoint',phase='update',last_source_row=cursor,merge=self.get('merges')+1)
        with self.db:
            size=self.db.execute('SELECT count(*) FROM vocab').fetchone()[0]
            self.db.execute('INSERT OR IGNORE INTO vocab VALUES(?,?,?,?)',
                            (size,selected['smiles'],count_atom(selected['smiles']),selected['frequency']))
            self.set('stage','choose')
            self.set('cursor',0)
            self.set('merges',self.get('merges')+1)
        self.notify('merge',merges=self.get('merges'),vocab_size=self.db.execute('SELECT count(*) FROM vocab').fetchone()[0],
                    selected=selected['smiles'],frequency=selected['frequency'])

    def train(self,engine,smiles):
        started=time.monotonic()
        self.spool(smiles)
        with mp.Pool(self.workers,initializer=_init_worker) as pool:
            while self.get('stage')!='complete':
                stage=self.get('stage')
                if stage=='initialize': self.initialize_graphs(pool)
                elif stage=='choose': self.choose()
                elif stage=='remove': self.remove_counts()
                elif stage=='update': self.update_graphs(pool)
                else: raise TokenizerError(f'Unknown disk checkpoint stage: {stage}')
        # Python's stable size sort preserves the legacy initial/learned order.
        vocab=list(self.db.execute('SELECT smiles,atoms,frequency FROM vocab ORDER BY atoms DESC,position'))
        engine.vocab_node=[s for s,_,_ in vocab]
        engine.vocab_node_index_map={s:i for i,s in enumerate(engine.vocab_node)}
        engine.vocab_node_stats={s:[a,f] for s,a,f in vocab}
        engine.vocab_node_indexed=[get_indexed_smiles(s,True) for s in engine.vocab_node]
        engine.initial_rings=self.get('initial_rings')
        engine.max_node_type=len(vocab)
        engine.max_atom_in_token=max(a for _,a,_ in vocab)
        self.notify('node_training_complete',seconds=time.monotonic()-started,vocab_size=len(vocab))

    def accepted_samples(self):
        # A separate read connection lets the multiprocessing feeder consume it.
        # The training checkpoint is stable by this point.
        with sqlite3.connect(self.directory/'npe.sqlite3') as reader:
            for (smiles,) in reader.execute('SELECT i.smiles FROM molecules m JOIN inputs i ON i.id=m.id ORDER BY m.id'):
                yield smiles

    def source_row(self,accepted_index):
        return self.db.execute('SELECT id FROM molecules WHERE accepted_index=?',(accepted_index,)).fetchone()[0]

    def rejections(self):
        yield from self.db.execute('SELECT i.id,i.smiles,r.reason FROM rejected r JOIN inputs i ON i.id=r.id ORDER BY i.id')
