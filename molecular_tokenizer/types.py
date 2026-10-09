"""Common model-facing payloads: no source SMILES are used for reconstruction."""
from dataclasses import asdict, dataclass
from typing import Literal


class TokenizerError(ValueError):
    pass


class UnsupportedMoleculeError(TokenizerError):
    pass


class FidelityError(UnsupportedMoleculeError):
    pass


class NotTrainedError(TokenizerError):
    pass


@dataclass(frozen=True)
class Connection:
    source: int
    target: int
    bond_type: int  # upstream convention: single=0, double=1, triple=2
    source_attachment: int
    target_attachment: int


@dataclass(frozen=True)
class MoleculeEncoding:
    backend: Literal['safe', 'npe']
    tokenizer_id: str
    token_ids: tuple[int, ...]  # sequence IDs for SAFE; motif node IDs for NPE
    connections: tuple[Connection, ...] = ()
    fidelity_checked: bool = False

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, value: dict) -> 'MoleculeEncoding':
        expected = {'backend', 'tokenizer_id', 'token_ids', 'connections', 'fidelity_checked'}
        if set(value) != expected:
            raise TokenizerError('Malformed encoding fields')
        result = cls(value['backend'], value['tokenizer_id'], tuple(value['token_ids']),
                     tuple(Connection(**edge) for edge in value['connections']), value['fidelity_checked'])
        result.validate()
        return result

    def validate(self) -> None:
        if self.backend not in {'safe', 'npe'} or not isinstance(self.tokenizer_id, str):
            raise TokenizerError('Invalid backend or tokenizer ID')
        if type(self.fidelity_checked) is not bool:
            raise TokenizerError('fidelity_checked must be boolean')
        if not self.token_ids or any(type(i) is not int or i < 0 for i in self.token_ids):
            raise TokenizerError('token_ids must be nonempty nonnegative integers')
        if self.backend == 'safe' and self.connections:
            raise TokenizerError('SAFE sequence encodings do not have graph connections')
        pairs = set()
        for edge in self.connections:
            if not isinstance(edge, Connection) or any(type(v) is not int for v in asdict(edge).values()):
                raise TokenizerError('Invalid connection fields')
            if not (0 <= edge.source < edge.target < len(self.token_ids)):
                raise TokenizerError('Connections require 0 <= source < target < node count')
            if edge.bond_type not in {0, 1, 2} or min(edge.source_attachment, edge.target_attachment) < 0:
                raise TokenizerError('Invalid bond type or attachment position')
            pair = (edge.source, edge.target)
            if pair in pairs:
                raise TokenizerError('Parallel motif connections are unsupported by this NPE adapter')
            pairs.add(pair)


@dataclass(frozen=True)
class TrainingReport:
    backend: str
    molecules: int
    vocab_size: int
    ring_vocab_size: int = 0
    edge_vocab_size: int = 0
    rejected_molecules: int = 0
    fragmentation: str = ''
    representation: str = ''
    motif_vocab_size: int = 0
    bpe_scope: str = ''
