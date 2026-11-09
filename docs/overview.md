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

[Machine-readable results, pseudocode and example IDs](vocabulary-results.json)
