"""Independent fragmentation and representation choices, with legacy codecs."""
from collections import Counter
from dataclasses import replace
import json

from rdkit import Chem
from tokenizers import Regex
from tokenizers.pre_tokenizers import Split
from safe import SAFEConverter

from .backends import SAFEBackend, NPEBackend, canonical_molecule
from .fragmentation import NPEPartitioner, brics_processor
from .types import TokenizerError, UnsupportedMoleculeError, TrainingReport
from ._vendor.graphbpe_utils import get_indexed_smiles, smiles_to_molecule


def collect_samples(smiles, *, skip_invalid=False, on_reject=None, on_reject_with_index=None,
                    reject_dummy=False):
    samples, indices, rejected = [], [], 0
    for index, sample in enumerate(smiles, 1):
        try:
            canonical_molecule(sample)
            if reject_dummy and '*' in sample:
                raise UnsupportedMoleculeError('Open dummy attachment points are unsupported')
            mol = Chem.MolFromSmiles(sample)
            Chem.Kekulize(mol, clearAromaticFlags=True)
        except Exception as error:
            if not skip_invalid:
                raise
            rejected += 1
            reason = f'{type(error).__name__}: {error}'
            if on_reject_with_index:
                on_reject_with_index(sample, reason, index)
            elif on_reject:
                on_reject(sample, reason)
            continue
        samples.append(sample)
        indices.append(index)
    if not samples:
        raise TokenizerError('Training corpus is empty')
    return samples, indices, rejected


class FragmentSAFEBackend(SAFEBackend):
    """BRICS or learned NPE cuts, then SAFE notation and sequence BPE.

    SAFE owns original atom properties and can veto stereo-unsafe cuts.
    ``fragment`` scope allows BPE to merge across atoms but never across dots.
    ``lexical`` retains the original atom/syntax pre-tokenization.
    """
    def __init__(self, fragmentation='npe', bpe_scope='fragment'):
        if fragmentation not in {'brics', 'npe'}:
            raise TokenizerError('fragmentation must be brics or npe')
        if bpe_scope not in {'lexical', 'fragment'}:
            raise TokenizerError('bpe_scope must be lexical or fragment')
        super().__init__(slicer='brics')
        self.fragmentation = fragmentation
        self.bpe_scope = bpe_scope
        self.npe = NPEBackend() if fragmentation == 'npe' else None
        self._configure()

    def _configure(self):
        if self.bpe_scope == 'fragment':
            self.tokenizer.pre_tokenizer = Split(Regex(r'\.'), behavior='isolated')
        if self.npe is not None:
            self.converter = SAFEConverter(slicer=NPEPartitioner(self.npe.engine), ignore_stereo=False)

    @property
    def config(self):
        return {'fragmentation': self.fragmentation, 'bpe_scope': self.bpe_scope}

    @property
    def files(self):
        return ('sequence.json',) + (NPEBackend.files if self.npe else ())

    def train(self, smiles, vocab_size, *, motif_vocab_size=350, ring_vocab_size=300,
              npe_model=None, num_workers=1, skip_invalid=False, on_reject=None,
              on_reject_with_index=None, min_frequency=2):
        if type(num_workers) is not int or num_workers < 1:
            raise TokenizerError('num_workers must be >= 1')
        if vocab_size < 95 or type(min_frequency) is not int or min_frequency < 1:
            raise TokenizerError('SAFE requires vocab_size >= 95 and min_frequency >= 1')
        if self.npe is None:
            if npe_model is not None:
                raise TokenizerError('npe_model requires NPE fragmentation')
            report = super().train(smiles, vocab_size, num_workers=num_workers,
                                   skip_invalid=skip_invalid, on_reject=on_reject,
                                   on_reject_with_index=on_reject_with_index, min_frequency=min_frequency)
        else:
            samples, indices, rejected = collect_samples(
                smiles, skip_invalid=skip_invalid, on_reject=on_reject,
                on_reject_with_index=on_reject_with_index, reject_dummy=True)
            if npe_model is not None:
                # Validate the source manifest/checksums before copying its partition model.
                from .tokenizer import MolecularTokenizer
                source = MolecularTokenizer.load(npe_model)
                backend = source.backend
                if isinstance(backend, FragmentSAFEBackend) and backend.npe is not None:
                    self.npe = backend.npe
                elif isinstance(backend, NPEBackend) and not isinstance(backend, BRICSDemoDiffBackend):
                    self.npe = backend
                else:
                    raise TokenizerError('npe_model must contain learned NPE motifs')
            else:
                if (type(motif_vocab_size) is not int or motif_vocab_size < 1
                        or type(ring_vocab_size) is not int or not 0 <= ring_vocab_size <= motif_vocab_size):
                    raise TokenizerError('Require motif_vocab_size > 0 and 0 <= ring_vocab_size <= motif_vocab_size')
                # Partition training does not need a DemoDiff edge vocabulary or decoder.
                self.npe.engine.train_node(samples, vocab_len=motif_vocab_size,
                                           vocab_ring_len=ring_vocab_size, num_processors=num_workers)
                if self.npe.vocab_size > motif_vocab_size:
                    raise TokenizerError('motif_vocab_size is smaller than the initial element/ring vocabulary')
            self._configure()

            def reject(sample, reason, filtered_index):
                if on_reject_with_index:
                    on_reject_with_index(sample, reason, indices[filtered_index - 1])
                elif on_reject:
                    on_reject(sample, reason)

            # Sequence conversion is serial for the learned callback; num_workers
            # controls NPE graph training. Never build an untrained worker slicer.
            report = super().train(samples, vocab_size, num_workers=1,
                                   min_frequency=min_frequency, skip_invalid=skip_invalid,
                                   on_reject_with_index=reject)
            report = replace(report, rejected_molecules=report.rejected_molecules + rejected,
                             ring_vocab_size=len(self.npe.engine.initial_rings),
                             motif_vocab_size=self.npe.vocab_size)
        return replace(report, fragmentation=self.fragmentation, representation='safe',
                       bpe_scope=self.bpe_scope)

    def state(self):
        state = {'sequence': SAFEBackend.state(self)}
        if self.npe:
            state['partition_model'] = self.npe.state()
        return state

    def save(self, directory):
        (directory / 'sequence.json').write_text(json.dumps(SAFEBackend.state(self)), encoding='utf-8')
        if self.npe:
            self.npe.save(directory)

    def load(self, directory):
        SAFEBackend.load(self, directory)
        if self.npe:
            self.npe.load(directory)
        self._configure()


class BRICSDemoDiffBackend(NPEBackend):
    """Rule-based partition with DemoDiff node/attachment representation.

    A finite observed fragment vocabulary; unseen motifs/edges are rejected.
    Inherits the upstream graph codec's fidelity limitations, checked by facade.
    """
    def __init__(self, fragmentation='brics'):
        if fragmentation != 'brics':
            raise TokenizerError('Expected brics fragmentation')
        super().__init__()
        self.engine.partition_processor = brics_processor

    @property
    def config(self):
        return {'fragmentation': 'brics'}

    def train(self, smiles, vocab_size, *, skip_invalid=False, on_reject=None,
              on_reject_with_index=None, num_workers=1):
        if num_workers != 1:
            raise TokenizerError('BRICS/DemoDiff vocabulary collection currently requires num_workers=1')
        samples, _, rejected = collect_samples(smiles, skip_invalid=skip_invalid,
            on_reject=on_reject, on_reject_with_index=on_reject_with_index)
        frequencies = Counter()
        for sample in samples:
            smi = canonical_molecule(sample)
            mol = smiles_to_molecule(smi, True)
            frequencies.update(brics_processor(smi, mol).subgraphs_smiles.values())
        if len(frequencies) > vocab_size:
            raise TokenizerError(f'BRICS requires {len(frequencies)} observed motifs; increase vocab_size (no silent truncation)')
        e = self.engine
        e.vocab_node = sorted(frequencies, key=lambda s: (-frequencies[s], s))
        e.vocab_node_index_map = {s:i for i,s in enumerate(e.vocab_node)}
        e.vocab_node_indexed = [get_indexed_smiles(s, True) for s in e.vocab_node]
        e.vocab_node_stats = {s:[smiles_to_molecule(s, True).GetNumAtoms(), frequencies[s]] for s in e.vocab_node}
        e.max_node_type = len(e.vocab_node)
        e.max_atom_in_token = max(v[0] for v in e.vocab_node_stats.values())
        for sample in samples:
            e.encode(sample, update_vocab_edge=True)
        e.unknown.clear()
        return TrainingReport('npe', len(samples), self.vocab_size,
                              edge_vocab_size=len(e.vocab_edge), rejected_molecules=rejected,
                              fragmentation='brics', representation='demodiff', motif_vocab_size=self.vocab_size)
