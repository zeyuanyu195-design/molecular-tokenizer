"""Separate representation algorithms behind a shared lifecycle."""
from abc import ABC, abstractmethod
import json
import multiprocessing as mp
from pathlib import Path
from typing import Iterable

from rdkit import Chem
from .types import Connection, TokenizerError, UnsupportedMoleculeError, TrainingReport


_SAFE_WORKER = None


def _init_safe_worker(slicer: str):
    global _SAFE_WORKER
    from rdkit import RDLogger
    RDLogger.DisableLog('rdApp.*')
    _SAFE_WORKER = SAFEBackend(slicer=slicer)


def _safe_worker_representation(smiles: str):
    try:
        return smiles, _SAFE_WORKER.representation(smiles), None
    except Exception as error:
        return smiles, None, f'{type(error).__name__}: {error}'


def canonical_molecule(smiles: str) -> str:
    if not isinstance(smiles, str) or not smiles or smiles != smiles.strip():
        raise TokenizerError('Expected a nonempty SMILES string without outer whitespace')
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise TokenizerError(f'Invalid SMILES: {smiles!r}')
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)


class TokenizerBackend(ABC):
    """Train may update vocabulary; encode/decode use a frozen vocabulary."""
    name: str
    files: tuple[str, ...]

    @property
    @abstractmethod
    def config(self) -> dict: ...

    @property
    @abstractmethod
    def vocab_size(self) -> int: ...

    @abstractmethod
    def train(self, smiles: Iterable[str], vocab_size: int, **options) -> TrainingReport: ...

    @abstractmethod
    def encode(self, smiles: str) -> tuple[tuple[int, ...], tuple[Connection, ...]]: ...

    @abstractmethod
    def decode(self, ids: tuple[int, ...], edges: tuple[Connection, ...]) -> str: ...

    @abstractmethod
    def state(self) -> dict: ...

    @abstractmethod
    def save(self, directory: Path) -> None: ...

    @abstractmethod
    def load(self, directory: Path) -> None: ...


class SAFEBackend(TokenizerBackend):
    """Official SAFE notation + SAFE lexical splitting + Rust BPE.

    This is a tokenization backend, not the pretrained SAFE-GPT model tokenizer.
    """
    name = 'safe'
    files = ('sequence.json',)

    def __init__(self, slicer: str = 'brics'):
        from safe import SAFEConverter, SAFESplitter
        from tokenizers import Tokenizer, decoders
        from tokenizers.models import BPE
        from tokenizers.pre_tokenizers import PreTokenizer
        if not isinstance(slicer, str) or slicer not in SAFEConverter.SUPPORTED_SLICERS:
            raise TokenizerError('Use a named SAFE slicer supported by SAFEConverter')
        self.slicer = slicer
        self.converter = SAFEConverter(slicer=slicer, ignore_stereo=False)
        # SAFE permits RDKit's extended 3–5 digit ring closures; the upstream
        # default splitter only recognizes 2-digit closures.
        splitter_pattern = SAFESplitter.REGEX_PATTERN.replace(r'\%[0-9]{2}', r'\%[0-9]{2,5}')
        self.splitter = SAFESplitter(pattern=splitter_pattern)
        self.tokenizer = Tokenizer(BPE(unk_token='[UNK]'))
        self.tokenizer.pre_tokenizer = PreTokenizer.custom(self.splitter)
        self.tokenizer.decoder = decoders.Fuse()

    @property
    def config(self):
        return {'slicer': self.slicer}

    @property
    def vocab_size(self):
        return self.tokenizer.get_vocab_size()

    def representation(self, smiles: str) -> str:
        if not isinstance(smiles, str) or not smiles or smiles != smiles.strip():
            raise TokenizerError('Expected a nonempty SMILES string without outer whitespace')
        # In SMILES, a wildcard/dummy atom is written with '*'. Avoid parsing the
        # molecule here: SAFEConverter already parses and validates it internally.
        if '*' in smiles:
            raise UnsupportedMoleculeError('Open dummy attachment points require a dedicated SAFE fragment workflow')
        # Ring-free and unsliceable molecules remain valid single SAFE blocks.
        try:
            return self.converter.encoder(smiles, canonical=True, randomize=False, allow_empty=True)
        except Exception as error:
            raise UnsupportedMoleculeError(f'SAFE conversion failed: {error}') from error

    def train(self, smiles, vocab_size, *, min_frequency=2, skip_invalid=False, on_reject=None,
              on_reject_with_index=None, num_workers=1):
        from tokenizers.trainers import BpeTrainer
        # Complete notation alphabet prevents silent [UNK] loss on unseen isotope numbers etc.
        if vocab_size < 95 or type(min_frequency) is not int or min_frequency < 1:
            raise TokenizerError('SAFE requires vocab_size >= 95 and min_frequency >= 1')
        if type(num_workers) is not int or num_workers < 1:
            raise TokenizerError('num_workers must be >= 1')
        count = rejected = 0

        def serial_representations():
            nonlocal count, rejected
            for sample_index, sample in enumerate(smiles, 1):
                try:
                    value = self.representation(sample)
                except Exception as error:
                    if not skip_invalid:
                        raise
                    rejected += 1
                    if on_reject_with_index is not None:
                        on_reject_with_index(sample, f'{type(error).__name__}: {error}', sample_index)
                    elif on_reject is not None:
                        on_reject(sample, f'{type(error).__name__}: {error}')
                    continue
                count += 1
                yield value

        def parallel_representations(pool):
            nonlocal count, rejected
            for sample_index, (sample, value, error) in enumerate(
                    pool.imap(_safe_worker_representation, smiles, chunksize=32), 1):
                if error is not None:
                    if not skip_invalid:
                        raise TokenizerError(error)
                    rejected += 1
                    if on_reject_with_index is not None:
                        on_reject_with_index(sample, error, sample_index)
                    elif on_reject is not None:
                        on_reject(sample, error)
                    continue
                count += 1
                yield value

        trainer = BpeTrainer(vocab_size=vocab_size, min_frequency=min_frequency,
                             special_tokens=['[UNK]'], initial_alphabet=[chr(i) for i in range(33, 127)],
                             show_progress=False)
        if num_workers == 1:
            self.tokenizer.train_from_iterator(serial_representations(), trainer=trainer)
        else:
            # SAFE conversion is CPU-bound and independent per SMILES. Keep the
            # Rust BPE trainer in the parent process and stream ordered worker
            # results so training remains memory-bounded and rejection records
            # preserve their source row association.
            with mp.Pool(num_workers, initializer=_init_safe_worker, initargs=(self.slicer,)) as pool:
                self.tokenizer.train_from_iterator(parallel_representations(pool), trainer=trainer)
        if count == 0:
            raise TokenizerError('Training corpus is empty')
        return TrainingReport(self.name, count, self.vocab_size, rejected_molecules=rejected)

    def encode(self, smiles):
        result = self.tokenizer.encode(self.representation(smiles), add_special_tokens=False)
        if self.tokenizer.token_to_id('[UNK]') in result.ids:
            raise UnsupportedMoleculeError('SAFE BPE produced an unknown token')
        return tuple(result.ids), ()

    def decode(self, ids, edges):
        if edges or any(i >= self.vocab_size for i in ids):
            raise TokenizerError('Invalid SAFE IDs or connections')
        if self.tokenizer.token_to_id('[UNK]') in ids:
            raise UnsupportedMoleculeError('Cannot decode unknown tokens')
        representation = self.tokenizer.decode(list(ids), skip_special_tokens=False)
        # Model outputs must decode as supplied; no automatic molecule repair/dummy deletion.
        result = self.converter.decoder(representation, canonical=True, fix=False, remove_dummies=False)
        return canonical_molecule(result)

    def state(self):
        # Custom Python pretokenizers are not serializable; restore them explicitly on load.
        pretok = self.tokenizer.pre_tokenizer
        try:
            self.tokenizer.pre_tokenizer = None
            return json.loads(self.tokenizer.to_str())
        finally:
            self.tokenizer.pre_tokenizer = pretok

    def save(self, directory):
        (directory / 'sequence.json').write_text(json.dumps(self.state(), ensure_ascii=False), encoding='utf-8')

    def load(self, directory):
        from tokenizers import Tokenizer
        from tokenizers.pre_tokenizers import PreTokenizer
        self.tokenizer = Tokenizer.from_file(str(directory / 'sequence.json'))
        self.tokenizer.pre_tokenizer = PreTokenizer.custom(self.splitter)


class NPEBackend(TokenizerBackend):
    """DemoDiff GraphBPE with frozen edge vocabulary and explicit graph connectivity."""
    name = 'npe'
    files = ('graph.node', 'graph.edge', 'graph.ring')

    def __init__(self):
        from ._vendor.graphbpe import MolecularGraphTokenizer
        self.engine = MolecularGraphTokenizer(kekulize=True, name='generic-npe')

    @property
    def config(self):
        return {}

    @property
    def vocab_size(self):
        return len(self.engine.vocab_node)

    def train(self, smiles, vocab_size, *, ring_vocab_size=300, num_workers=1,
              skip_invalid=False, on_reject=None, on_reject_with_index=None):
        if type(num_workers) is not int or num_workers < 1:
            raise TokenizerError('num_workers must be >= 1')
        if type(ring_vocab_size) is not int or not 0 <= ring_vocab_size <= vocab_size:
            raise TokenizerError('ring_vocab_size must be between 0 and vocab_size')
        samples = []
        rejected = 0
        for sample_index, sample in enumerate(smiles, 1):  # Upstream NPE revisits the molecule set for iterative merging.
            try:
                canonical_molecule(sample)
                samples.append(sample)
            except Exception as error:
                if not skip_invalid:
                    raise
                rejected += 1
                if on_reject_with_index is not None:
                    on_reject_with_index(sample, f'{type(error).__name__}: {error}', sample_index)
                elif on_reject is not None:
                    on_reject(sample, f'{type(error).__name__}: {error}')
        if not samples:
            raise TokenizerError('Training corpus is empty')
        # The GraphBPE parser validates each sample again as it builds a graph;
        # avoid a redundant full corpus-wide RDKit pass here.
        self.engine.train_node(samples, vocab_len=vocab_size, vocab_ring_len=ring_vocab_size,
                               num_processors=num_workers)
        if self.vocab_size > vocab_size:
            raise TokenizerError('vocab_size is smaller than the initial atom/ring vocabulary')
        for sample in samples:
            nodes, _ = self.engine.encode(sample, update_vocab_edge=True)
            if any(i >= self.vocab_size for i in nodes):
                raise UnsupportedMoleculeError('Training molecule has a motif absent from the learned vocabulary')
        self.engine.unknown.clear()
        return TrainingReport(self.name, len(samples), self.vocab_size,
                              len(self.engine.initial_rings), len(self.engine.vocab_edge), rejected)

    def encode(self, smiles):
        canonical_molecule(smiles)
        self.engine.unknown.clear()
        try:
            nodes, adj, edge_existence = self.engine.encode(
                smiles, update_vocab_edge=False, return_edge_existence=True)
            if any(i >= self.vocab_size for i in nodes):
                raise UnsupportedMoleculeError('Molecule has an unknown NPE motif or atom')
            connections = []
            for i in range(len(nodes)):
                for j in range(i + 1, len(nodes)):
                    a, b = adj[i][j], adj[j][i]
                    if a < 0 and b < 0:
                        if edge_existence[i][j]:
                            raise UnsupportedMoleculeError('Attachment type is absent from the frozen edge vocabulary')
                        continue
                    if a < 0 or b < 0:
                        raise UnsupportedMoleculeError('Attachment type is absent from the frozen edge vocabulary')
                    left, right = self.engine.vocab_edge[a], self.engine.vocab_edge[b]
                    if len(left) != 2 or len(right) != 2:
                        raise UnsupportedMoleculeError('Parallel bonds between motifs are unsupported by upstream decode')
                    if left[1] != right[1]:
                        raise TokenizerError('Inconsistent directed bond types')
                    connections.append(Connection(i, j, left[1], left[0], right[0]))
            return tuple(nodes), tuple(connections)
        finally:
            self.engine.unknown.clear()

    def decode(self, ids, edges):
        if any(i >= self.vocab_size for i in ids):
            raise UnsupportedMoleculeError('Unknown NPE node token')
        size = len(ids)
        bond_adj = [[-1] * size for _ in ids]
        positions = [[-1] * size for _ in ids]
        from ._vendor.graphbpe_utils import smiles_to_molecule
        atom_counts = [smiles_to_molecule(self.engine.vocab_node_indexed[i], True).GetNumAtoms() for i in ids]
        for edge in edges:
            if (edge.source_attachment >= atom_counts[edge.source]
                    or edge.target_attachment >= atom_counts[edge.target]):
                raise TokenizerError('Attachment position is outside the motif')
            if ((edge.source_attachment, edge.bond_type) not in self.engine.vocab_edge
                    or (edge.target_attachment, edge.bond_type) not in self.engine.vocab_edge):
                raise UnsupportedMoleculeError('Connection is absent from the frozen edge vocabulary')
            i, j = edge.source, edge.target
            bond_adj[i][j] = bond_adj[j][i] = edge.bond_type
            positions[i][j], positions[j][i] = edge.source_attachment, edge.target_attachment
        value, _, unknown = self.engine.decode(list(ids), bond_adj, positions, replace_unknown_with_random=False)
        if value is None or unknown:
            raise UnsupportedMoleculeError('NPE graph could not be decoded into a valid molecule')
        return canonical_molecule(value)

    def state(self):
        e = self.engine
        return {'node': e.vocab_node, 'indexed': e.vocab_node_indexed,
                'stats': e.vocab_node_stats, 'edges': e.vocab_edge,
                'edge_stats': [[list(key), list(e.vocab_edge_stats[key])] for key in e.vocab_edge],
                'rings': e.initial_rings, 'kekulize': e.kekulize}

    def save(self, directory):
        self.engine.save(str(directory / 'graph'))
        # Upstream omits .edge when no inter-motif bonds were learned.
        (directory / 'graph.edge').touch(exist_ok=True)

    def load_prefix(self, prefix):
        prefix = Path(prefix)
        self.engine.load(str(prefix))
        # Empty ring/edge files are valid for an atom-only or wholly merged corpus.
        if not Path(str(prefix) + '.ring').read_text(encoding='utf-8').strip():
            self.engine.initial_rings = []
        edge_path = Path(str(prefix) + '.edge')
        if not edge_path.exists() or not edge_path.read_text(encoding='utf-8').strip():
            self.engine.vocab_edge = []
            self.engine.vocab_edge_stats = {}
            self.engine.max_edge_type = 0
        if not self.engine.kekulize:
            raise TokenizerError('This NPE adapter requires a kekulized vocabulary')

    def load(self, directory):
        self.load_prefix(directory / 'graph')
