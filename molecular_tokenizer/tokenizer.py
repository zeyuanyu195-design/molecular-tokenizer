"""Generic facade with transactional training and vocabulary-bound encodings."""
from dataclasses import replace
import hashlib
import json
import os
from pathlib import Path
import shutil
import tempfile
from threading import RLock
from typing import Iterable, Iterator

from .backends import SAFEBackend, NPEBackend, TokenizerBackend, canonical_molecule
from .types import FidelityError, MoleculeEncoding, NotTrainedError, TokenizerError, TrainingReport

BACKENDS = {'safe': SAFEBackend, 'npe': NPEBackend}


def build_backend(name, config):
    if 'fragmentation' in config:
        from .pipelines import FragmentSAFEBackend, BRICSDemoDiffBackend
        return (FragmentSAFEBackend if name == 'safe' else BRICSDemoDiffBackend)(**config)
    return BACKENDS[name](**config)


class MolecularTokenizer:
    """Shared train/encode/decode/save/load API; different payload topology by backend.

    strict=True rejects encodings that change the canonical isomeric molecule.
    An encoding contains only learned IDs/connectivity, never the original SMILES.
    """
    def __init__(self, backend: str | None = None, *, strict: bool = True,
                 fragmentation: str | None = None, representation: str | None = None,
                 bpe_scope: str | None = None, **backend_config):
        if backend is None:
            fragmentation = fragmentation or 'brics'
            representation = representation or 'safe'
            if fragmentation not in {'brics', 'npe'} or representation not in {'safe', 'demodiff'}:
                raise TokenizerError('Choose fragmentation=brics/npe and representation=safe/demodiff')
            if backend_config:
                raise TokenizerError('Use explicit fragmentation/representation options without legacy backend options')
            if representation == 'safe':
                backend = 'safe'
                backend_config = {'fragmentation': fragmentation, 'bpe_scope': bpe_scope or 'fragment'}
            else:
                if bpe_scope is not None:
                    raise TokenizerError('bpe_scope is only meaningful for SAFE representation')
                backend = 'npe'
                backend_config = {'fragmentation': 'brics'} if fragmentation == 'brics' else {}
        elif representation is not None or bpe_scope is not None or fragmentation is not None:
            raise TokenizerError('Do not mix legacy backend with fragmentation/representation/bpe_scope')
        if backend not in BACKENDS:
            raise TokenizerError(f'Unknown backend {backend!r}; use safe or npe')
        if type(strict) is not bool:
            raise TokenizerError('strict must be boolean')
        self.backend: TokenizerBackend = build_backend(backend, backend_config)
        self.strict = strict
        self._trained = False
        self._tokenizer_id = None
        self._lock = RLock()

    @property
    def vocab_size(self):
        return self.backend.vocab_size

    @property
    def fragmentation(self):
        return self.backend.config.get('fragmentation', self.backend.config.get('slicer', 'npe'))

    @property
    def representation(self):
        return 'safe' if self.backend.name == 'safe' else 'demodiff'

    @property
    def tokenizer_id(self):
        self._require_trained()
        return self._tokenizer_id

    def _require_trained(self):
        if not self._trained:
            raise NotTrainedError('Train or load a vocabulary before encoding/decoding')

    def _fingerprint(self, backend=None):
        backend = self.backend if backend is None else backend
        state = {'schema': 1, 'backend': backend.name, 'config': backend.config,
                 'vocabulary': backend.state()}
        return hashlib.sha256(json.dumps(state, sort_keys=True, separators=(',', ':')).encode()).hexdigest()

    def train(self, smiles: Iterable[str], vocab_size: int, **options) -> TrainingReport:
        if isinstance(smiles, (str, bytes)) or type(vocab_size) is not int or vocab_size < 1:
            raise TokenizerError('Provide an iterable of SMILES and a positive integer vocab_size')
        with self._lock:
            candidate = build_backend(self.backend.name, self.backend.config)
            report = candidate.train(smiles, vocab_size, **options)
            fingerprint = self._fingerprint(candidate)
            # Preserve an existing model if training raises; swap only after a full successful fit.
            self.backend = candidate
            self._trained = True
            self._tokenizer_id = fingerprint
            return replace(report,
                           fragmentation=report.fragmentation or self.fragmentation,
                           representation=report.representation or self.representation,
                           motif_vocab_size=report.motif_vocab_size or
                           (self.vocab_size if self.representation == 'demodiff' else 0))

    def to_safe(self, smiles: str) -> str:
        """Return the SAFE string, checking molecular identity in strict mode."""
        with self._lock:
            self._require_trained()
            if self.representation != 'safe':
                raise TokenizerError('to_safe requires SAFE representation')
            original = canonical_molecule(smiles)
            value = self.backend.representation(smiles)
            if self.strict:
                decoded = self.backend.converter.decoder(value, canonical=True, fix=False, remove_dummies=False)
                if canonical_molecule(decoded) != original:
                    raise FidelityError('SAFE conversion changes the molecule')
            return value

    def encode(self, smiles: str) -> MoleculeEncoding:
        with self._lock:
            self._require_trained()
            original = canonical_molecule(smiles)
            ids, edges = self.backend.encode(smiles)
            result = MoleculeEncoding(self.backend.name, self.tokenizer_id, ids, edges)
            result.validate()
            if self.strict:
                reconstructed = self.backend.decode(ids, edges)
                if canonical_molecule(reconstructed) != original:
                    raise FidelityError(f'{self.backend.name} changes molecule: {smiles!r} -> {reconstructed!r}')
                result = replace(result, fidelity_checked=True)
            return result

    def encode_many(self, smiles: Iterable[str]) -> Iterator[MoleculeEncoding]:
        for sample in smiles:
            yield self.encode(sample)

    def decode(self, encoding: MoleculeEncoding) -> str:
        with self._lock:
            self._require_trained()
            if not isinstance(encoding, MoleculeEncoding):
                raise TokenizerError('decode expects MoleculeEncoding, including graph connectivity for NPE')
            encoding.validate()
            if encoding.backend != self.backend.name or encoding.tokenizer_id != self.tokenizer_id:
                raise TokenizerError('Encoding belongs to a different backend or vocabulary')
            return self.backend.decode(encoding.token_ids, encoding.connections)

    def save(self, directory: str | Path) -> Path:
        with self._lock:
            self._require_trained()
            target = Path(directory).resolve()
            if target.exists():
                raise FileExistsError(f'Choose a new model directory: {target}')
            target.parent.mkdir(parents=True, exist_ok=True)
            temporary = Path(tempfile.mkdtemp(prefix='tokenizer-', dir=target.parent))
            try:
                self.backend.save(temporary)
                checksums = {name: hashlib.sha256((temporary / name).read_bytes()).hexdigest()
                             for name in self.backend.files}
                manifest = {'schema_version': 1, 'backend': self.backend.name,
                            'backend_config': self.backend.config, 'strict': self.strict,
                            'tokenizer_id': self.tokenizer_id, 'sha256': checksums}
                (temporary / 'tokenizer.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
                os.replace(temporary, target)
            finally:
                if temporary.exists():
                    shutil.rmtree(temporary)
            return target

    @classmethod
    def load(cls, directory: str | Path) -> 'MolecularTokenizer':
        target = Path(directory)
        manifest = json.loads((target / 'tokenizer.json').read_text(encoding='utf-8'))
        if manifest.get('schema_version') != 1:
            raise TokenizerError('Unsupported tokenizer schema')
        model = cls(manifest['backend'], strict=manifest['strict'])
        model.backend = build_backend(manifest['backend'], manifest['backend_config'])
        if set(manifest['sha256']) != set(model.backend.files):
            raise TokenizerError('Model manifest files do not match the backend')
        for name in model.backend.files:
            if hashlib.sha256((target / name).read_bytes()).hexdigest() != manifest['sha256'][name]:
                raise TokenizerError(f'Model checksum mismatch: {name}')
        model.backend.load(target)
        model._trained = True
        model._tokenizer_id = model._fingerprint()
        if model.tokenizer_id != manifest['tokenizer_id']:
            raise TokenizerError('Loaded vocabulary does not match the saved fingerprint')
        return model

    @classmethod
    def from_graphbpe(cls, prefix: str | Path, *, strict=True) -> 'MolecularTokenizer':
        model = cls('npe', strict=strict)
        model.backend.load_prefix(prefix)
        model._trained = True
        model._tokenizer_id = model._fingerprint()
        return model
