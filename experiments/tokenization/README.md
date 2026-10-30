# Tokenization experiments

This directory archives the latest **original SMILES → original DemoDiff** experiment and the helper scripts used to investigate fixed SAFE boundaries. The root README and GitHub Pages remain the results presentation; execution details live here.

## Published results and counting units

| Experiment / cohort | Molecules | Total tokens | Mean |
|---|---:|---:|---:|
| Original SMILES → DemoDiff, successful records | 49,998 / 50,000 | 428,440 | 8.569143 |
| GraphBPE+SAFE, all exact-boundary records | 44,329 | 372,957 | 8.413386 |
| Original DemoDiff, matched cohort | 44,328 | 372,924 | 8.412832 |
| GraphBPE+SAFE fixed motifs, matched cohort | 44,328 | 372,924 | 8.412832 |
| BRICS+SAFE text tokenizer, matched cohort | 44,328 | 2,669,535 | 60.222320 |
| GraphBPE+SAFE text control, matched cohort | 44,328 | 3,051,888 | 68.847861 |

DemoDiff and fixed-motif lengths count **node/motif IDs only**, excluding edge arrays, attachment positions and SAFE connectivity syntax. SAFE text lengths include dots and connection labels, with no special tokens, padding or truncation. These are different units, not a direct comparison of complete lossless sequence lengths.

Of 49,992 GraphBPE+SAFE records, 5,663 have changed boundaries and are excluded from the exact-boundary cohort. Do not assume every change has the same stereochemical cause. One exact-boundary record, CHEMBL1213498, is absent from BRICS, leaving 44,328 matched molecules. The matched cohorts have zero per-molecule length disagreements. Encoding success alone is not a full round-trip validation of all 49,998 DemoDiff outputs.

- [Recomputed statistics, errors, cohort exclusions and source hashes](results/token_lengths.json)
- [Original DemoDiff run summary (portable paths)](results/demodiff_original_smiles.summary.json)
- [Three checked examples](results/examples.json)
- [Artifact sizes, SHA-256 hashes and publication scope](results/artifact_manifest.json)
- [Environment used for this verification](results/verification_environment.json)

## Scripts

| Script | Purpose |
|---|---|
| `export_demodiff_tokens.py` | Run upstream DemoDiff on original SMILES; save one row per input, including errors. This produced the latest baseline. |
| `export_exact_graphbpe.py` | Filter existing GraphBPE+SAFE records by validated exact boundaries and reuse recorded node IDs; does not segment again. |
| `demodiff_fixed_safe.py` | Map existing SAFE dot-fragments to compatible DemoDiff motif keys without calling GraphBPE encode/merge. Unknown keys are explicit. |
| `summarize_token_lengths.py` | Recompute full and matched cohorts from saved outputs, validating identities, partitions and counts. |
| `legacy/direct_demodiff_safe.py` | Earlier experiment: complete SAFE strings → original DemoDiff, which may resegment. **Not the original-SMILES baseline.** |
| `experimental/safe_fragment_tokenizer.py` | Alternative fragment dictionary. Literal mode preserves full fragment strings; rings mode additionally needs ring-label metadata for decoding. This is dictionary construction, not frequency-based BPE. |
| `experimental/compare_fixed_safe_bpe.py` | Alternative shared text BPE trained on both SAFE representations, with merges constrained by dots. Does not enforce one motif = one ID. |

The `experimental/` approaches did **not** produce the published table. The mapper's IDs are in SAFE fragment order; recorded GraphBPE IDs are in `bpe_groups` order. Node IDs alone do not preserve attachment positions or connection labels.

## Environment and upstream dependency

Run from the repository root with Python 3.11 or newer. DemoDiff export and mapping need `rdkit` and `tqdm`; shared text BPE needs `tokenizers`; statistics additionally need `numpy` and this project's SAFE tokenizer dependencies.

```sh
python -m pip install -e ".[model]"
git clone https://github.com/liugangcode/DemoDiff.git ../DemoDiff
git -C ../DemoDiff checkout 0acb85ef446c175a3e3a60bbad43de97d97ca76f
```

Keep the upstream checkout separate. The scripts import its `downstream/graphbpe.py` and `utils/token.py`, and load `data/tokenizer/vocab3000ring300/pretrain-token.{node,edge,ring}`. They do not replace upstream code or retrain the vocabulary. SHA-256 hashes in the run summary identify the actual source/vocabulary bytes. Recorded environment versions describe local verification, not a tested cross-platform lockfile.

## Reproduce the latest baseline

Replace paths to match your machine. The input is the **original ChEMBL TSV**, column `canonical_smiles`; do not substitute `safe.txt`. `--limit 50000` means the first 50,000 data rows, excluding the header. Without `--limit`, the exporter reads the whole input.

```sh
python experiments/tokenization/export_demodiff_tokens.py --demodiff-root ../DemoDiff --input safe/data/chembl_37_chemreps.txt --smiles-column canonical_smiles --limit 50000 --output runs/demodiff_original_smiles/smiles_50000_tokens.jsonl
```

The output contains `record_index`, `chembl_id`, `smiles`, `status`, node IDs/counts and edge/bond/attachment arrays. Error rows are retained and excluded from the mean. The companion `.summary.json` records counts, source and vocabulary hashes. Existing outputs are refused.

Core algorithm:

```text
load the frozen upstream DemoDiff vocabulary
for each of the first 50,000 original SMILES rows:
    try:
        node_ids, adjacency = upstream.encode(SMILES, update_vocab_edge=False)
        bonds, positions = upstream.get_bond_position_by_vocab(adjacency)
        save identifiers, SMILES, node_ids, adjacency, bonds, positions
        count len(node_ids)
    except:
        save identifiers, SMILES and error
mean = sum(successful node counts) / number of successful records
```

## Fixed-boundary helpers and matched statistics

These commands reuse the existing GraphBPE+SAFE and BRICS+SAFE records. The exact-boundary exporter requires `bpe_boundary_exact`, `bpe_groups`, `safe_groups`, `bpe_node_ids`, fragment counts and `roundtrip_ok` from the original encoding pipeline; plain SAFE strings alone are insufficient to establish equality with the original atom partition.

```sh
python experiments/tokenization/export_exact_graphbpe.py --input safe/data/output/records.jsonl --output-dir runs/graphbpe_exact
python experiments/tokenization/demodiff_fixed_safe.py --demodiff-root ../DemoDiff --safe "c12nc3oc1C.N34CC5.C5#N.C2.N4"
python experiments/tokenization/summarize_token_lengths.py --graph safe/data/output/records.jsonl --brics safe/data/output_brics/records.jsonl --demodiff runs/demodiff_original_smiles/smiles_50000_tokens.jsonl --text-tokenizer expts/tokenizer/tokenizer-custom.json --output runs/token_lengths.json
```

Use `--demodiff-root` for the supported upstream mapper example; the older `--graphbpe-dir` option expects a standalone module layout, not this repository's `safe` package.

```text
exact = records whose saved GraphBPE and SAFE atom partitions agree
validate fragment counts, saved IDs and recorded molecular round-trip flags
reuse saved GraphBPE node IDs; retain the SAFE string for connectivity
matched = exact ∩ BRICS records ∩ successful original-SMILES DemoDiff records
check record identifiers and original SMILES before comparing lengths
load the SAFE text tokenizer with its custom pre-tokenizer
report full-cohort and matched-cohort lengths separately
```

Validation command (small fixtures, not another 50,000-molecule run):

```sh
python -m unittest discover -s experiments/tokenization -p "test_*.py"
```

All helper CLIs support `--help`. The legacy and experimental scripts retain their separate outputs and should not be substituted into `summarize_token_lengths.py`.

## Large artifact publication

This code update includes small summaries and three examples. It does **not** add the full latest per-record outputs or the full input corpus. Existing repository data paths already use Git LFS; their availability depends on those LFS objects having been uploaded.

For a fixed experiment snapshot, publish compressed JSONL outputs as versioned GitHub Release assets and include the manifest, source commit, vocabulary hashes and commands. For data that must follow every checkout, use Git LFS. Do not add these large JSONL files directly to ordinary Git history. The manifest's `included_in_this_commit` describes this publication only; it does not claim that older LFS artifacts are absent from GitHub.
