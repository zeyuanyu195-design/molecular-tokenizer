# SAFE connection walkthrough

The existing results page now includes a six-step walkthrough: input graph,
recorded GraphBPE partition, accepted boundaries, SAFE connections, fixed-fragment
vocabulary lookup, and reconstruction. Previous/next, direct step selection, and
optional five-second playback all use the same checked-in trace.

| Saved experiment record | GraphBPE motifs | SAFE fragments | Accepted / proposed cuts | Result |
| --- | ---: | ---: | ---: | --- |
| CHEMBL13214, input record 10047 | 8 | 8 | 7 / 7 | Exact boundaries; all eight motif IDs found |
| CHEMBL55653, input record 50000 | 4 | 3 | 2 / 3 | C11=C12 retained with STEREOE; merged F1 is absent from the node vocabulary |

Atom indices are zero-based indices in the saved converter callback molecule.
The changed example is excluded from the fixed-boundary comparison, not from
the demonstration: it reconstructs successfully. Its three SAFE fragments must
not be assigned the original four GraphBPE IDs.

## What carries connectivity

For CHEMBL55653, the actual string is
`c12ncc3cc1/C=C/c1ccccc1.N=C(N)N2.Cl3`.
The pair `2` links the aromatic carbon at source atom 4 to nitrogen at atom 3;
the pair `3` links aromatic carbon 7 to chlorine 8. Both are single bonds.
The `1` pairs close internal rings and are not cross-fragment connections.
The `/C=C/` region remains within F1, preserving the specified double-bond
geometry. The slashes are directional single-bond syntax adjacent to the double
bond. The actual double-bond stereo was checked as STEREOE.

The converter creates matching dummy endpoints on each accepted cut, writes
fragments and replaces each matching dummy pair with a ring-closure label.
Attachment atoms are specified by the label positions in the string; bond and
stereochemical syntax are retained. A dot is a separator, not an edge list.
Selecting a label highlights both occurrences, their atoms and their bond.
The same label may be reused after a pair closes; the UI tracks occurrences.

The graph, colors and atom mapping in the JSON are **visualization metadata**.
They are not decoder inputs. `SAFEConverter.decoder(complete_safe_string, ...)`
was called without an attachment table. Motif IDs alone are not a reversible
representation in the current implementation.

## Reproduce and validate

From a checkout with the project's chemistry dependencies installed:

```bash
python experiments/tokenization/build_connection_demo.py \
  --repo-root . \
  --vocab-prefix /path/to/DemoDiff/data/tokenizer/vocab3000ring300/pretrain-token \
  --output docs/connection-demo.json
```

The input is the two genuine saved records in `connection_demo_records.json`.
The builder replays the current pipeline and fails if SAFE strings, GraphBPE
groups/IDs, proposed cuts or accepted cuts differ. It verifies:

- source and complete SAFE have equal canonical isomeric graph identities;
- a stereochemistry-preserving graph isomorphism aligns SAFE atoms with the
  recorded source partition (symmetric equivalent mappings can exist);
- every cross-fragment label pair is a real accepted cut and vice versa;
- the boundary-preserving example's lookup IDs equal its recorded motif IDs,
  after accounting for SAFE fragment order;
- the complete SAFE decoder output matches the original isomeric molecule.

The generated JSON includes source/vocabulary/input fingerprints, RDKit version,
paired token occurrences, fragment membership, coordinates, SVGs and validation
flags. Tested with RDKit 2025.09.1. The displayed boundary policy belongs to this
repository's converter. No intermediate GraphBPE merge animation is fabricated;
step 2 shows the final recorded partition. No full 50,000-record rerun is needed.
