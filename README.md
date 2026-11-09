# Reversible GraphBPE + SAFE tokenization

[Interactive results and token-by-token example](https://zeyuanyu195-design.github.io/molecular-tokenizer/#vocabulary-results)

## Latest result — 44,329 boundary-preserving records

We compare **whole-fragment lookup with base-token fallback** against **whole-fragment lookup with within-fragment BPE fallback**, at the same total model vocabulary budget. Both preserve complete SAFE strings, including connection labels, bond/stereo symbols and dots.

| Total vocabulary | Records | Whole-fragment mean | Hybrid mean | Length reduction | Whole-fragment total | Hybrid total |
|---:|---:|---:|---:|---:|---:|---:|
| 16,000 | 44,329 | 35.088542 | 19.919376 | 43.23% | 1,555,440 | 883,006 |
| 32,000 | 44,329 | 27.402220 | 18.183695 | 33.64% | 1,214,713 | 806,065 |

Reduction = 1 − hybrid total / whole-fragment total. No BOS/EOS, padding or truncation. These are complete encoding lengths, not the earlier motif-only counts (~8.41) that exclude connectivity.

**Scope:** all 44,329 records have identical GraphBPE and final SAFE partitions. The 5,663 changed-boundary records are excluded. The 44,329 records contain 44,328 distinct original canonical-SMILES strings; the table counts records without deduplication. This is different from the earlier 44,328-record cohort matched across BRICS and DemoDiff.

## From graph fragments to a bounded vocabulary

1. Start with the existing SMILES → GraphBPE partition → SAFE serialization. Do not rerun GraphBPE during text tokenization.
2. Keep the existing split: **35,518 training / 4,268 validation / 4,543 test records**. The same non-stereochemical molecular graph and its stereoisomers stay in one split.
3. Keep **M = 2,136** base entries: 1,880 existing SAFE entries + 256 byte-fallback entries. Count literal complete fragments on TRAIN only; e.g. `C4c1cccs1` and `C8c1cccs1` are different entries.
4. Whole-fragment scheme: add the most frequent K−M complete fragments. A miss uses base text tokens; dots remain separate tokens.
5. Hybrid scheme: allocate part of K−M to complete fragments, and the rest to TRAIN-only BPE substrings learned inside each fragment. Start from base tokens, never merge across dots, never rewrite labels. Identical payloads share an output ID.
6. Compare 25%, 50%, 75% whole-fragment allocations on VALIDATION, then freeze the selected configuration. A complete-fragment hit emits one ID; a miss uses within-fragment BPE.
7. Decode by concatenating global vocabulary byte payloads. The reconstructed complete SAFE carries connectivity; no original SAFE or per-molecule edge table is needed.

| Scheme / total budget | Base entries | Dedicated whole-fragment entries | Added substring entries |
|---|---:|---:|---:|
| Whole / 16k | 2,136 | 13,864 | 0 |
| Hybrid / 16k | 2,136 | 6,932 | 6,932 |
| Whole / 32k | 2,136 | 29,864 | 0 |
| Hybrid / 32k | 2,136 | 7,466 | 22,398 |

BPE can itself learn a complete-fragment string; the categories describe how entries were allocated, not exclusive chemical types. A substring need not be a valid standalone molecule.

## Pseudocode

### Whole-fragment scheme

```text
TRAIN(train_SAFE, K):
    base = existing_SAFE_vocabulary + 256_byte_fallback_entries
    M = 2136
    fragments = split_each_training_SAFE_on_dots()
    ranked = frequency_rank(fragments, lexical_tie_break)
    remove_strings_already_represented_by_one_base_token(ranked)
    add_top(K - M, ranked)  # literal labels retained
    return vocabulary

ENCODE(safe, vocabulary):
    for fragment in safe.split("."):
        if complete_fragment_entry_exists(fragment):
            emit its_ID
        else:
            emit base_tokenize(fragment)
            # Use UTF-8 byte fallback if exact reconstruction fails.
        emit dot_ID between fragments
    return IDs

DECODE(IDs, vocabulary):
    safe = UTF8_decode(concatenate vocabulary.byte_payload[ID])
    return safe
    # Parse the COMPLETE SAFE string to recover the molecule.
    # No original SAFE side table or per-molecule edge table.
```

### Hybrid scheme

```text
TRAIN(train_SAFE, K, whole_fraction):
    keep the same base vocabulary M = 2136
    reserve floor((K-M) * whole_fraction) frequent whole fragments
    base_encode each TRAIN fragment separately
    learn BPE merges weighted by TRAIN fragment frequency
    # Never merge across dots; preserve connection labels.
    append merge-prefix substrings until total model IDs == K
    reuse an existing ID for an identical byte payload
    return vocabulary, within_fragment_BPE

SELECT_ON_VALIDATION():
    compare whole fractions 25%, 50%, 75% of (K-M)
    minimize validation total length; freeze before test evaluation

ENCODE(safe, vocabulary, within_fragment_BPE):
    for fragment in safe.split("."):
        if complete_fragment_entry_exists(fragment):
            emit its_ID
        else:
            base_IDs = lossless_base_encode(fragment)
            emit within_fragment_BPE(base_IDs).model_IDs
        emit dot_ID between fragments
    return IDs

DECODE(IDs, vocabulary):
    safe = UTF8_decode(concatenate vocabulary.byte_payload[ID])
    return safe
    # Parse the COMPLETE SAFE string to recover the molecule.
    # No original SAFE side table or per-molecule edge table.
```

## CHEMBL45183 — 39 → 23 tokens at 16k

Original record_index: 39224 (TEST). All eight SAFE fragment boundaries are preserved.

```text
C123CCC(C1)C45C26.C7(=O)C8NC(=O)[C@H]9N.C8c1cccs1.C9C(=O)O.O67.C3.C4.C5
```

| Fragment | Old tokens | Hybrid tokens | Hybrid text payloads (one quoted string per token) |
|---|---:|---:|---|
| `C123CCC(C1)C45C26` | 1 | 5 | `C123` · `CCC(C` · `1)` · `C45` · `C26` |
| `C7(=O)C8NC(=O)[C@H]9N` | 17 | 4 | `C7(=O)` · `C8` · `NC(=O)[C@H]9` · `N` |
| `C8c1cccs1` | 9 | 2 | `C8` · `c1cccs1` |
| `C9C(=O)O` | 1 | 1 | `C9C(=O)O` |
| `O67` | 1 | 1 | `O67` |
| `C3` | 1 | 1 | `C3` |
| `C4` | 1 | 1 | `C4` |
| `C5` | 1 | 1 | `C5` |
| Separating dots | 7 | 7 | Each `.` has ID 11 |
| **Total** | **39** | **23** | |

Fragment 1 increases from 1 to 5 tokens because it loses its dedicated whole-fragment slot; fragments 2 and 3 decrease from 17 to 4 and 9 to 2. Net saving: 16 tokens. `C8c1cccs1` now uses `C8` and `c1cccs1`, preserving both external label 8 and internal ring label 1.

Final sequence (line breaks only separate original fragments):

```text
"C123", "CCC(C", "1)", "C45", "C26", ".",
"C7(=O)", "C8", "NC(=O)[C@H]9", "N", ".",
"C8", "c1cccs1", ".",
"C9C(=O)O", ".",
"O67", ".",
"C3", ".",
"C4", ".",
"C5"
```

## Validation and interpretation

- Both budgets and both methods restored **44,329 / 44,329** primary SAFE strings exactly in the prior per-record checks. The new report filters and recomputes those verified records; it does not retrain the vocabulary.
- At 16k: 32,854 records become shorter, 4,597 unchanged, 6,878 longer. At 32k: 23,127 shorter, 13,442 unchanged, 7,760 longer.
- The primary table includes TRAIN and VALIDATION records; it is a descriptive full-cohort result, not a held-out test score. Only new entries were fit on TRAIN; provenance of the pre-existing base/NPE vocabularies has not been audited.
- On the 4,543-record TEST set, whole → hybrid means are **38.781862 → 20.742461 (16k)** and **33.608849 → 19.193044 (32k)**. The test set was explored in earlier diagnosis, so it is not a pristine final holdout.
- Pure within-fragment BPE controls achieve TEST means **20.744662 (16k)** and **19.193044 (32k)**. Most demonstrated compression comes from substring reuse; dedicated whole-fragment slots have not shown a substantial extra compression advantage.
- Existing fragment boundaries are preserved, but one motif can produce multiple model tokens. Lossless encoding does not establish valid diffusion generation or attachment-label pairing in newly generated molecules.

[Machine-readable results, pseudocode and example IDs](docs/vocabulary-results.json)

<details>
<summary>Historical experiments: BRICS, original DemoDiff and motif-only counts</summary>

The material below describes earlier experiments only; its implementation limitations and sample counts do not describe the new reversible codecs above.

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

GraphBPE segmentation, SAFE serialization, molecular reconstruction, and fixed-boundary motif lookup are implemented for this comparison. GraphBPE+SAFE reconstruction here retains the complete SAFE string; DemoDiff reconstruction uses separate edge and attachment-position information.

The three displayed examples have verified round trips. The full DemoDiff batch has not undergone a complete molecular round-trip validation. Equal motif counts do not by themselves establish equal information content, full encoding compression, or lossless generation.

## Evidence

- [Machine-readable results and the three examples](docs/results.json)
- [GraphBPE / SAFE boundary recording](safe/batch_chembl.py)
- [Current SAFE stereochemistry handling](safe/converter.py)
- [SAFE text-token counting](safe/compare_safe_token_lengths.py)
- [Original DemoDiff implementation](https://github.com/liugangcode/DemoDiff/blob/0acb85ef446c175a3e3a60bbad43de97d97ca76f/downstream/graphbpe.py)

Source-result and tokenizer fingerprints are included in the machine-readable results. Statistics describe the local experimental outputs, not a new training or evaluation run performed by this webpage.

</details>
