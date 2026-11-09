# Whole-fragment and hybrid token walkthrough

The main page's existing graph-to-SAFE walkthrough now has eight steps:
input graph, recorded GraphBPE motifs, boundary validation, paired SAFE labels,
complete-fragment dictionary lookup, fallback encoding, actual ID sequences,
and ID-to-payload reconstruction. The molecular drawing and selectable connection
labels remain available. No unrecorded GraphBPE merge history is simulated.

The three saved records in `codec_demo_records.json` are genuine experiment
inputs. Two preserve boundaries; the third is an excluded stress example.

| Example | Motifs → SAFE fragments | Whole → hybrid, 16k | Whole → hybrid, 32k |
| --- | ---: | ---: | ---: |
| CHEMBL45183 / record 39224 | 8 → 8 | 39 → 23 | 39 → 18 |
| CHEMBL13214 / record 10047 | 8 → 8 | 26 → 18 | 15 → 17 |
| CHEMBL55653 / record 50000 | 4 → 3 | 27 → 7 | 27 → 7 |

These are complete text-token lengths, including dots and connection syntax.
Compression is not guaranteed for every example. The changed-boundary example
does not contribute to the 44,329-record comparison.

## What lookup and fallback mean

Lookup uses the entire literal SAFE fragment, including labels. An added whole
entry or an exact base entry emits a single ID and skips fallback. Otherwise,
the whole scheme uses the saved base tokenizer (byte fallback is available),
while hybrid applies its frozen BPE within that fragment's base-token sequence.
The two dictionaries allocate different numbers of whole-fragment entries at
the same total budget, so a hit in one is not necessarily a hit in the other.

For CHEMBL45183 at 16k, `C123CCC(C1)C45C26` has whole-scheme ID 14286.
It is absent from the hybrid complete-fragment dictionary. Hybrid encodes it as
`C123 | CCC(C | 1) | C45 | C26`, IDs
`11140 | 10047 | 9447 | 2191 | 12524`.
In particular, hybrid ID 9447 really has UTF-8 payload `3129` (hex), or `1)`.
It is a text substring, not an independent chemical motif. The same numeric ID
in a different scheme or budget need not have the same payload.

## Reproduction

Use the original frozen experiment artifacts, not freshly retrained vocabularies.
From the repository root, with chemistry and tokenizer dependencies installed:

```bash
python experiments/tokenization/build_connection_demo.py \
  --repo-root . \
  --records experiments/tokenization/codec_demo_records.json \
  --vocab-prefix /path/to/DemoDiff/data/tokenizer/vocab3000ring300/pretrain-token \
  --output /path/to/graphs.json

python experiments/tokenization/build_codec_walkthrough.py \
  --graphs /path/to/graphs.json \
  --whole-dir /path/to/bounded_safe_vocab_20261106 \
  --hybrid-dir /path/to/hybrid_safe_vocab_20261106 \
  --output docs/connection-demo.json
```

The graph builder checks replayed SAFE strings, partitions, cuts, atom mappings
and label pairs. The codec builder reloads each frozen encoder and decoder and
checks all 12 example/budget/scheme combinations:

- Every displayed ID resolves in the actual dictionary; whole hits stay single.
- No fallback token crosses a dot; concatenated fragment payloads equal the input.
- Full ID sequences agree with the production codec calls.
- Decoding the IDs produces exactly the original complete SAFE string.
- The real SAFE decoder receives that restored text alone, with `fix=False`;
  its canonical isomeric molecule equals the original.

The JSON includes full-vocabulary SHA-256 fingerprints and an excerpt of the
shared dictionary entries used by these examples. The browser validates all
ID sequences again and reconstructs text incrementally from IDs and dictionary
entries. The original SAFE is only an equality target, not a decoding input.
Graph metadata is only for visualization. Molecular reconstruction is checked
offline; the browser does not execute RDKit or SAFE. Valid round trips do not
establish the validity of arbitrary sequences sampled by a generative model.
