# Molecular tokenization: GraphBPE, SAFE and DemoDiff

[Explore the interactive results](https://zeyuanyu195-design.github.io/molecular-tokenizer/)

## Main finding

On **44,328 matched molecules**, fixed-boundary GraphBPE+SAFE and the original DemoDiff tokenizer produce **identical motif/node counts for every molecule**: **372,924 tokens in total**, averaging **8.412832 per molecule**. DemoDiff receives the original `canonical_smiles`, not a SAFE string.

These are motif/node counts **excluding connectivity overhead**. They are not directly comparable as complete encoding costs to SAFE text-token lengths, which include connectivity symbols.

## Train, encode and decode

| Method | Train / vocabulary used | Encode | Decode | Counted unit |
| --- | --- | --- | --- | --- |
| BRICS+SAFE | Fixed BRICS rules; existing SAFE text tokenizer | SMILES → BRICS → SAFE → text IDs | Text IDs → SAFE → SMILES | SAFE text tokens |
| Original DemoDiff | NPE/GraphBPE motif learning; existing node/ring/edge vocabulary | Original SMILES → GraphBPE → node IDs + edge/attachment information | Node IDs + bond adjacency + attachment positions → SMILES | Node IDs only |
| GraphBPE+SAFE | Same existing GraphBPE node/ring vocabulary | SMILES → GraphBPE → SAFE → fixed-boundary motif lookup | Complete SAFE → SMILES | Motif IDs only |

No vocabulary is retrained in this experiment. The graph vocabulary contains 3,000 node entries and 300 ring entries; its original training corpus and that of the SAFE text tokenizer have not been verified. The additional frequency-based vocabulary previously proposed is **not** part of these results.

## Token-length results

### Shared molecules

All rows below use the same **44,328 molecules**, matched by record index and checked against ChEMBL ID and original SMILES.

| Method | Molecules | Total tokens | Mean | Median | P95 | Maximum |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| BRICS+SAFE | 44,328 | 2,669,535 | 60.222320 | 52 | 107 | 1,487 |
| Original DemoDiff | 44,328 | 372,924 | 8.412832 | 5 | 22 | 223 |
| GraphBPE+SAFE · fixed motifs | 44,328 | 372,924 | 8.412832 | 5 | 22 | 223 |
| GraphBPE+SAFE · text control | 44,328 | 3,051,888 | 68.847861 | 56 | 133 | 2,779 |

BRICS+SAFE and the GraphBPE+SAFE text control use the **same SAFE text tokenizer**, with no special tokens, padding or truncation, and with dropout disabled. Only these two rows measure the same text-token unit. The two motif rows exclude connectivity information.

### Full available cohorts

These cohorts have different inclusion criteria and are not matched-sample comparisons. GraphBPE+SAFE fixed motifs includes only the 44,329 boundary-preserving records.

| Method | Molecules | Total tokens | Mean | Median | P95 | Maximum |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| BRICS+SAFE | 49,997 | 3,024,042 | 60.484469 | 52 | 109 | 1,487 |
| Original DemoDiff | 49,998 | 428,440 | 8.569143 | 5 | 23 | 246 |
| GraphBPE+SAFE · fixed motifs | 44,329 | 372,957 | 8.413386 | 5 | 22 | 223 |
| GraphBPE+SAFE · text control | 49,992 | 3,435,594 | 68.722876 | 56 | 135 | 2,779 |

The original DemoDiff experiment processes the first **50,000 original SMILES**: **49,998 succeed**, **2 fail SMILES parsing**, and **0 successfully encoded rows report unknown nodes**.

## Boundary preservation and sample selection

- GraphBPE+SAFE succeeds for 49,992 of the 50,000 input records.
- 5,663 successful records have final SAFE boundaries different from the original GraphBPE partition and are excluded from the fixed-boundary comparison.
- 44,329 retain the original boundaries. Their motif counts also match original-SMILES DemoDiff counts record by record: 372,957 total, mean 8.413386.
- One of those records, CHEMBL1213498, lacks a successful BRICS result, leaving 44,328 shared records for the three-method table.

The current SAFE implementation restricts some cuts to preserve stereochemistry. Boundary changes are not limited to cutting stereogenic double bonds: omitted cuts also involve single, aromatic and triple bonds. For CHEMBL55653, an omitted double-bond cut changes four original GraphBPE motifs into three SAFE fragments while preserving the reconstructed molecule.

SAFE dots separate fragments; **paired labels, their attachment positions and bond symbols**, rather than dots alone, encode cross-fragment connectivity.

## Core pseudocode

The training blocks describe the conceptual roles of `train`. The reported measurements load existing vocabularies. Encoding and decoding blocks distinguish the three evaluated paths.

### BRICS+SAFE

```text
TRAIN(training_SAFE_strings, vocabulary_size):
    use predefined BRICS fragmentation rules
    learn a SAFE text-token vocabulary from training strings
    return text_tokenizer

# This experiment loads an existing text tokenizer.
# Its original training corpus has not been verified.
```

```text
ENCODE(smiles, text_tokenizer):
    safe = SAFE.encode(smiles, slicer=BRICS)
    ids = text_tokenizer.encode(safe,
        add_special_tokens=False, padding=False,
        truncation=False, dropout=0)
    return safe, ids, len(ids)
```

```text
DECODE(text_ids, text_tokenizer):
    safe = text_tokenizer.decode(text_ids)
    return SAFE.decode(safe)
```

### Original DemoDiff

```text
TRAIN(training_molecular_graphs, target_size):
    initialize graph vocabulary and ring handling
    repeat until target_size or no eligible merge:
        count eligible adjacent-subgraph pairs
        merge the most frequent pair in training graphs
        add the resulting motif to the node vocabulary
    retain node, ring and edge/attachment vocabularies

# Conceptual NPE training, not rerun for this experiment.
# The original vocabulary training corpus is unverified.
```

```text
ENCODE(original_smiles, tokenizer):
    node_ids, adj = tokenizer.encode(
        original_smiles, update_vocab_edge=False)
    bonds, positions = tokenizer.get_bond_position_by_vocab(adj)
    return node_ids, bonds, positions, len(node_ids)

# Input is the original canonical_smiles, not a SAFE string.
```

```text
DECODE(node_ids, bonds, positions, tokenizer):
    return tokenizer.decode(node_ids, bonds, positions,
        replace_unknown_with_random=False)
```

### GraphBPE+SAFE

```text
TRAIN(training_molecular_graphs, target_size):
    return NPE_TRAIN(training_molecular_graphs, target_size)

# For the reported experiment, load the existing vocabulary.
# No additional SAFE-fragment vocabulary is trained.
```

```text
ENCODE(smiles, vocabulary):
    groups = GraphBPE.partition(smiles, vocabulary)
    safe = SAFE.serialize(smiles, proposed_partition=groups)
    if final_SAFE_partition(safe) != groups:
        exclude from the boundary-preserving comparison
    ids = []
    for fragment in SAFE_order(safe):
        subgraph = extract_from_complete_SAFE_graph(fragment)
        key = DemoDiff_compatible_motif_key(subgraph)
        ids.append(vocabulary.lookup(key))
    return safe, ids, len(ids)

# No second GraphBPE merge pass during fixed-boundary lookup.
# Aggregate lengths reuse len(record.bpe_node_ids) after
# exact-boundary validation; examples also verify this mapping.
```

```text
DECODE(safe, motif_ids):
    return SAFE.decode(safe)

# Connectivity is retained in the complete SAFE string.
# Decoding from motif_ids alone is not implemented.
```

## Three boundary-preserving examples

All three examples belong to the shared cohort and pass molecular reconstruction checks for all three methods. Fixed motif IDs are listed in **SAFE fragment order**; DemoDiff node IDs follow its own graph order. The same motif may appear more than once.

### 1. CHEMBL45973 — 5 GraphBPE motifs

Original SMILES:

```text
Cc1nc(N(N)CCC#N)oc1C
```

| Representation | Fragments / nodes | Counted token length |
| --- | ---: | ---: |
| BRICS+SAFE | 1 | 20 text tokens |
| Original DemoDiff | 5 | 5 node IDs |
| GraphBPE+SAFE | 5 | 5 motif IDs |
| GraphBPE+SAFE text control | 5 | 28 text tokens |

**Complete BRICS SAFE**

```text
Cc1nc(N(N)CCC#N)oc1C
```

**Complete GraphBPE SAFE**

```text
c12nc3oc1C.N34CC5.C5#N.C2.N4
```

<details>
<summary>SAFE fragments, motif keys and IDs</summary>

| SAFE fragment | DemoDiff-compatible motif key | ID |
| --- | --- | ---: |
| `c12nc3oc1C` | `CC1=CN=CO1` | 2577 |
| `N34CC5` | `CCN` | 2895 |
| `C5#N` | `C#N` | 2942 |
| `C2` | `C` | 2971 |
| `N4` | `N` | 2972 |

DemoDiff node IDs: `[2971, 2972, 2895, 2942, 2577]`.

</details>

### 2. CHEMBL13214 — 8 GraphBPE motifs

Original SMILES:

```text
CN(C)Cc1ccc(CSCCNC(=O)NCCSCc2ccc(CN(C)C)o2)o1
```

| Representation | Fragments / nodes | Counted token length |
| --- | ---: | ---: |
| BRICS+SAFE | 15 | 79 text tokens |
| Original DemoDiff | 8 | 8 node IDs |
| GraphBPE+SAFE | 8 | 8 motif IDs |
| GraphBPE+SAFE text control | 8 | 60 text tokens |

**Complete BRICS SAFE**

```text
c15ccc6o1.c17ccc8o1.CN%10C.C%15C%11.C39=O.C%12C%16.N%13(C)C.C%105.C%146.S%14%15.N3%11.N9%12.S%164.C47.C%138
```

**Complete GraphBPE SAFE**

```text
C3c1ccc4o1.C6(=O)NCC7.c18ccc9o1.C5CN6.C9N(C)C.CN3C.C4S5.S7C8
```

<details>
<summary>SAFE fragments, motif keys and IDs</summary>

| SAFE fragment | DemoDiff-compatible motif key | ID |
| --- | --- | ---: |
| `C3c1ccc4o1` | `CC1=CC=CO1` | 2509 |
| `C6(=O)NCC7` | `CCNC=O` | 2722 |
| `c18ccc9o1` | `C1=COC=C1` | 2682 |
| `C5CN6` | `CCN` | 2895 |
| `C9N(C)C` | `CN(C)C` | 2855 |
| `CN3C` | `CNC` | 2905 |
| `C4S5` | `CS` | 2945 |
| `S7C8` | `CS` | 2945 |

DemoDiff node IDs: `[2682, 2895, 2722, 2905, 2945, 2945, 2509, 2855]`.

</details>

### 3. CHEMBL34237 — 12 GraphBPE motifs

Original SMILES:

```text
CCCCNCC1OC(OC2C(N)CC(N)C(OC3OC(CN)C(O)C(O)C3N)C2O)C(O)C(N)C1O
```

| Representation | Fragments / nodes | Counted token length |
| --- | ---: | ---: |
| BRICS+SAFE | 9 | 79 text tokens |
| Original DemoDiff | 12 | 12 node IDs |
| GraphBPE+SAFE | 12 | 12 motif IDs |
| GraphBPE+SAFE text control | 12 | 84 text tokens |

**Complete BRICS SAFE**

```text
C1%10OC4C(O)C(N)C1O.C16C(N)CC(N)C7C1O.C15OC%11C(O)C(O)C1N.CCCC8.N89.C9%10.O46.O57.C%11N
```

**Complete GraphBPE SAFE**

```text
C14OC7C5C(N)C16.C18C(N)CC9C%10C1%14.C1%11OC%12C(O)C(O)C1%13.CCCCNC4.O78.O%10%11.C%12N.N9.N%13.O%14.O5.O6
```

<details>
<summary>SAFE fragments, motif keys and IDs</summary>

| SAFE fragment | DemoDiff-compatible motif key | ID |
| --- | --- | ---: |
| `C14OC7C5C(N)C16` | `NC1CCOCC1` | 2341 |
| `C18C(N)CC9C%10C1%14` | `NC1CCCCC1` | 2233 |
| `C1%11OC%12C(O)C(O)C1%13` | `OC1CCOCC1O` | 2136 |
| `CCCCNC4` | `CCCCNC` | 2607 |
| `O78` | `O` | 2973 |
| `O%10%11` | `O` | 2973 |
| `C%12N` | `CN` | 2940 |
| `N9` | `N` | 2972 |
| `N%13` | `N` | 2972 |
| `O%14` | `O` | 2973 |
| `O5` | `O` | 2973 |
| `O6` | `O` | 2973 |

DemoDiff node IDs: `[2973, 2972, 2973, 2972, 2973, 2973, 2973, 2940, 2233, 2341, 2607, 2136]`.

</details>

## Interpretation and implementation status

GraphBPE segmentation, SAFE serialization, molecular reconstruction, and fixed-boundary motif lookup are implemented for this comparison. **A reversible model-token sequence decoded from motif IDs alone has not been implemented.** GraphBPE+SAFE reconstruction here retains the complete SAFE string; DemoDiff reconstruction uses separate edge and attachment-position information.

The three displayed examples have verified round trips. The full DemoDiff batch has not undergone a complete molecular round-trip validation. Equal motif counts do not by themselves establish equal information content, full encoding compression, or lossless generation.

## Evidence

- [Machine-readable results and the three examples](results.json)
- [GraphBPE / SAFE boundary recording](../safe/batch_chembl.py)
- [Current SAFE stereochemistry handling](../safe/converter.py)
- [SAFE text-token counting](../safe/compare_safe_token_lengths.py)
- [Original DemoDiff implementation](https://github.com/liugangcode/DemoDiff/blob/0acb85ef446c175a3e3a60bbad43de97d97ca76f/downstream/graphbpe.py)

Source-result and tokenizer fingerprints are included in the machine-readable results. Statistics describe the local experimental outputs, not a new training or evaluation run performed by this webpage.
