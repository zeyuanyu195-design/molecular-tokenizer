# GraphBPE 与 BRICS：ChEMBL 批量切分及 token 长度比较

本说明针对 `E:\GraphBpe+SAFE\SAFE-main` 当前项目，使用已有的 Conda 环境 `graphbpe-safe`。运行顺序为：准备环境 → BRICS 切分 → GraphBPE 切分 → 比较 SAFE 文本的平均 token 长度。

本文命令使用 Windows CMD / Anaconda Prompt。Python 命令均为单行，终端自动折行不影响执行；不要添加 PowerShell 的反引号。所有数据和脚本参数使用绝对路径，因此无需先 `cd` 到项目目录。

## 1. 脚本与数据位置

| 用途 | 完整路径 |
| --- | --- |
| BRICS 批处理脚本 | `E:\GraphBpe+SAFE\SAFE-main\safe\batch_chembl_brics.py` |
| GraphBPE 批处理脚本 | `E:\GraphBpe+SAFE\SAFE-main\safe\batch_chembl.py` |
| token 长度比较脚本 | `E:\GraphBpe+SAFE\SAFE-main\safe\compare_safe_token_lengths.py` |
| ChEMBL 输入数据 | `E:\GraphBpe+SAFE\SAFE-main\safe\data\chembl_37_chemreps.txt` |
| GraphBPE 切分词表目录 | `E:\GraphBpe+SAFE\SAFE-main\safe\data\tokenizer\vocab3000ring300` |
| SAFE 文本 tokenizer | `E:\GraphBpe+SAFE\SAFE-main\expts\tokenizer\tokenizer-custom.json` |
| 依赖列表 | `E:\GraphBpe+SAFE\SAFE-main\requirements.txt` |

当前 `.txt` 输入文件实际上是带表头的 TSV（制表符分隔）。脚本至少要求 `chembl_id` 和 `canonical_smiles` 两列；当前文件还包含 `standard_inchi` 和 `standard_inchi_key`。脚本根据表头而非扩展名识别分隔符。

## 2. 激活环境

在 Anaconda Prompt 中执行：

```cmd
conda activate graphbpe-safe
python -c "import sys; print(sys.executable)"
```

当前机器应显示：

```text
E:\programming\conda\envs\graphbpe-safe\python.exe
```

如果需要按项目依赖文件安装或补齐依赖，执行：

```cmd
python -m pip install -r "E:\GraphBpe+SAFE\SAFE-main\requirements.txt"
```

已配置好的环境无需每次运行前重复安装。当前依赖文件内容为：

```txt
safe-mol==0.2.0
rdkit==2025.9.1
datamol==0.13.0
numpy==2.5.2
networkx==3.6.1
loguru==0.7.3
tqdm==4.70.1
torch==2.14.1
tokenizers==0.23.3
transformers==5.20.1
huggingface-hub==1.33.0
fsspec==2026.9.0
```

## 3. 运行 BRICS 切分

先处理输入文件前 1000 条数据记录：

```cmd
python -B "E:\GraphBpe+SAFE\SAFE-main\safe\batch_chembl_brics.py" --input "E:\GraphBpe+SAFE\SAFE-main\safe\data\chembl_37_chemreps.txt" --output-dir "E:\GraphBpe+SAFE\SAFE-main\safe\data\output_brics" --workers 1 --batch-size 256 --limit 1000 --delimiter auto --allow-stereo-coarsening
```

该脚本使用 RDKit BRICS 规则提出切分边，再由 SAFE 转换器生成 SAFE 字符串，并检查切分分区和分子重建结果。

输出目录会自动创建：

```text
E:\GraphBpe+SAFE\SAFE-main\safe\data\output_brics
```

## 4. 运行 GraphBPE 切分

使用同一输入文件的前 1000 条数据记录：

```cmd
python -B "E:\GraphBpe+SAFE\SAFE-main\safe\batch_chembl.py" --input "E:\GraphBpe+SAFE\SAFE-main\safe\data\chembl_37_chemreps.txt" --output-dir "E:\GraphBpe+SAFE\SAFE-main\safe\data\output" --vocab-prefix "E:\GraphBpe+SAFE\SAFE-main\safe\data\tokenizer\vocab3000ring300\pretrain-token" --workers 1 --batch-size 256 --limit 1000 --delimiter auto --allow-stereo-coarsening
```

该脚本加载已有的 在前5万条数据中训练的GraphBPE 词表，由 BPE 子图边界提出切分边，再调用 SAFE 转换器生成 SAFE 字符串。运行此命令不会训练新的 GraphBPE 词表。

`--vocab-prefix` 是三个文件共有的路径前缀，不加扩展名。当前所需文件为：

```text
E:\GraphBpe+SAFE\SAFE-main\safe\data\tokenizer\vocab3000ring300\pretrain-token.node
E:\GraphBpe+SAFE\SAFE-main\safe\data\tokenizer\vocab3000ring300\pretrain-token.edge
E:\GraphBpe+SAFE\SAFE-main\safe\data\tokenizer\vocab3000ring300\pretrain-token.ring
```

当前 `batch_chembl.py` 默认指向项目根目录下的 `data\tokenizer`，而实际词表位于 `safe\data\tokenizer`，因此应保留命令中的显式 `--vocab-prefix`。

输出目录会自动创建：

```text
E:\GraphBpe+SAFE\SAFE-main\safe\data\output
```

## 5. 批处理参数说明

| 参数 | 适用脚本 | 默认值 / 是否必填 | 含义 |
| --- | --- | --- | --- |
| `--input` | 两者 | 必填 | 带表头的 CSV/TSV 数据文件路径。至少包含 `chembl_id`、`canonical_smiles`。 |
| `--output-dir` | 两者 | 必填 | 结果保存目录；不存在时自动创建。 |
| `--vocab-prefix` | 仅 GraphBPE | 默认指向项目根目录的 `data\tokenizer\vocab3000ring300\pretrain-token` | GraphBPE 词表公共前缀；代码在后面追加 `.node`、`.edge`、`.ring`。本文命令显式指定实际路径。 |
| `--workers` | 两者 | `1` | 转换工作进程数，必须大于等于 1。`1` 表示串行处理；例如 `4` 使用 4 个工作进程。各工作进程持有自己的转换器，多进程会增加内存占用。 |
| `--batch-size` | 两者 | `256` | 主进程每次读取并提交处理的输入记录数，必须大于等于 1。它是数据处理批量大小，不是模型训练 batch size。 |
| `--limit` | 两者 | `0` | 最多读取多少条数据记录，不含表头。`1000` 为前 1000 条；`0` 为处理全部。失败记录也占用这个限额。 |
| `--delimiter` | 两者 | `auto` | `auto` 根据表头判断：含制表符则按 TSV，否则按逗号分隔。`tab` 强制 TSV，`comma` 强制 CSV。 |
| `--allow-stereo-coarsening` | 两者 | 默认关闭 | 无需填写值的开关。允许 SAFE 的立体保护调整原切分边界；仍要求最终分子重建及分区检查通过。 |
| `-h` / `--help` | 两者 | 可选 | 打印参数列表并退出，不执行切分。 |

命令中的 `-B` 是 Python 参数，用于禁止写入 `.pyc` 字节码缓存，与切分结果无关。

### 立体保护开关如何选择

本文两条命令均加入 `--allow-stereo-coarsening`，与当前已生成结果的设置一致。

- 不加该参数：要求算法提出的切分边与 SAFE 实际使用的切分边完全一致。若 SAFE 因立体保护改变切分边，该记录写入 `failed.jsonl`。
- 加上该参数：允许这种边界变化，可能使最终 SAFE 片段比原切分更粗；仍检查重建后的分子是否保持立体、荷电和同位素信息。该参数不等于忽略立体化学，也不保证每条记录都成功。
- 比较两种方法时，两条命令应采用相同的开关设置。若研究严格保留算法边界的结果，应从两条命令中同时删除该参数。

### 如何处理全量数据

将两条命令中的 `--limit 1000` 都改成 `--limit 0` 即可。需要并行时，可把 `--workers 1` 改成适合机器资源的进程数，例如 `--workers 4`。

同一个输出目录会被重复使用，`safe.txt`、`records.jsonl`、`failed.jsonl` 会在运行时重新写入，`summary.json` 在正常完成后重新写入。保留试验结果时，请给新运行设置不同的 `--output-dir`；中断运行后不要将目录中旧的 `summary.json` 当成本次完成标志。

## 6. 切分结果保存在哪里

按本文命令运行后，两组结果分别保存在：

```text
E:\GraphBpe+SAFE\SAFE-main\safe\data\output\
    safe.txt
    records.jsonl
    failed.jsonl
    summary.json

E:\GraphBpe+SAFE\SAFE-main\safe\data\output_brics\
    safe.txt
    records.jsonl
    failed.jsonl
    summary.json
```

| 文件 | 内容及用途 |
| --- | --- |
| `safe.txt` | 每行一条成功生成的 SAFE 字符串，可作为语料，也供后续 token 长度比较使用。失败记录不会写入。 |
| `records.jsonl` | 每行一个成功记录的 JSON 对象，包含原始记录编号、ChEMBL ID、原始 SMILES、SAFE 编码、切分边、片段分组及检查结果。 |
| `failed.jsonl` | 每行一个失败记录，包含原始标识、`error_type` 和 `error`。没有失败时文件可为空。 |
| `summary.json` | 整次运行的统计，包括处理/成功/失败数量、边界完全一致数量、错误类型、耗时、输入路径和版本。GraphBPE 还保存词表哈希。某些零计数字段可能不出现。 |

`safe.txt` 第 n 行对应同一目录 `records.jsonl` 第 n 条成功记录；如果存在失败，它不一定对应原始输入的第 n 条记录。追溯输入时使用 `record_index`（从 1 开始，不含表头）和 `chembl_id`。

成功记录中的主要切分字段：

| 字段 | 含义 |
| --- | --- |
| `safe` | 最终 SAFE 编码字符串。 |
| `requested_cuts` / `used_cuts` | 原算法提出的切分边 / SAFE 实际接受的切分边。 |
| `safe_groups` / `num_safe_fragments` | 实际 SAFE 片段的原子分组 / 片段数量。 |
| `bpe_groups` / `num_bpe_fragments` | GraphBPE 原始子图分组 / 子图数量，仅 GraphBPE 输出。 |
| `bpe_node_ids` | GraphBPE 子图的词表 ID，仅 GraphBPE 输出；不等于 SAFE 文本 tokenizer 的 ID 序列。 |
| `brics_groups` / `num_brics_fragments` | 原始 BRICS 切分分组 / 片段数量，仅 BRICS 输出。 |
| `bpe_boundary_exact` / `brics_boundary_exact` | 相应算法的提出切分边与实际切分边是否完全一致。 |
| `roundtrip_ok` | 是否通过分子重建一致性检查；写入成功记录时为 `true`。 |

## 7. 比较平均 token 长度

两种方法的切分都结束后，直接读取生成的两份 `safe.txt`：

```cmd
python -B "E:\GraphBpe+SAFE\SAFE-main\safe\compare_safe_token_lengths.py" --project-root "E:\GraphBpe+SAFE\SAFE-main" --tokenizer "E:\GraphBpe+SAFE\SAFE-main\expts\tokenizer\tokenizer-custom.json" --graphbpe "E:\GraphBpe+SAFE\SAFE-main\safe\data\output\safe.txt" --brics "E:\GraphBpe+SAFE\SAFE-main\safe\data\output_brics\safe.txt"
```

当前比较脚本已包含自动保存逻辑：会创建 `comparison` 目录、覆盖写入下列文件，并在终端打印结果和保存位置。无需手动添加 `>`，也没有 `--output` 参数。

```text
E:\GraphBpe+SAFE\SAFE-main\safe\data\comparison\results.json
```

若修改了切分命令的输出目录，也要相应修改比较命令中的 `--graphbpe` 和 `--brics`。比较结果仍保存到 `--project-root` 下的 `safe\data\comparison\results.json`。

| 比较参数 | 默认值 / 是否必填 | 含义 |
| --- | --- | --- |
| `--project-root` | 默认 `E:\GraphBpe+SAFE\SAFE-main` | 本地项目根目录，用于导入项目代码并确定比较结果保存位置。 |
| `--tokenizer` | 必填 | 两份语料共同使用的已训练 tokenizer JSON 文件，也可指定包含 `tokenizer.json` 的目录。 |
| `--graphbpe` | 必填 | GraphBPE 输出的 `safe.txt` 路径。 |
| `--brics` | 必填 | BRICS 输出的 `safe.txt` 路径。 |
| `--add-special-tokens` | 默认关闭 | 无需填写值的开关；添加后计入 tokenizer 配置的特殊 token，例如 `[CLS]`、`[SEP]`。 |
| `-h` / `--help` | 可选 | 查看参数帮助并退出。 |

统计逻辑：对每条非空 SAFE 文本使用同一个 tokenizer 编码，以 ID 列表长度为 token 数，最后计算 `mean_tokens = total_tokens / num_sequences`。默认不计特殊 token、不 padding、不截断，并关闭 BPE dropout。空行被计数后跳过；编码失败会报错终止，不会悄悄丢弃。

GraphBPE 的 `.node/.edge/.ring` 用于分子切分，`tokenizer-custom.json` 用于 SAFE 文本的 tokenization，两类词表不能互相替代。

## 8. 如何读比较结果

| JSON 字段 | 含义 |
| --- | --- |
| `graphbpe.mean_tokens` / `brics.mean_tokens` | 两种 SAFE 表示的平均 token 长度。 |
| `num_sequences` / `total_tokens` | 纳入统计的非空序列数量 / token 总数。 |
| `min_tokens` / `max_tokens` | 单条序列 token 数的最小值 / 最大值。 |
| `unk_tokens` / `sequences_with_unk` | 未知 token 总数 / 含未知 token 的序列数量。 |
| `comparison.graphbpe_minus_brics` | GraphBPE 均值减去 BRICS 均值；负值表示 GraphBPE 更短。 |
| `comparison.graphbpe_reduction_percent` | `(BRICS均值 - GraphBPE均值) / BRICS均值 × 100`；正值表示 GraphBPE 更短，负值表示更长。 |
| `comparison.same_sequence_count` | 两份文件的非空序列数是否相同，不证明分子集合相同。 |
| `comparison.same_molecule_set_verified` | 当前脚本固定写为 `false`，因为它没有读取分子标识进行匹配。 |
| `tokenizer_sha256` | 文本 tokenizer 文件的哈希，用于标识所用词表文件。 |

两种切分方法可能有不同的失败记录。正式比较应使用共同成功的同一批分子；针对同一个输入文件，可按两份 `records.jsonl` 的 `record_index` 对齐并核对分子标识，再统计对应的 `safe` 字段。当前比较脚本不会自动求交集。

平均 token 长度也依赖 tokenizer 的训练语料与词表。报告结果时应说明 tokenizer 来源及特殊 token 设置，不能将某个文本 tokenizer 下的长度直接等同于 GraphBPE 子图数量。

### 已有结果快照

2026-10-19 读取当前 `comparison\results.json` 得到以下结果。这是已有文件的快照，不代表本文重新执行了切分；重新运行会更新结果。

| 指标 | GraphBPE | BRICS |
| --- | ---: | ---: |
| 序列数量 | 1000 | 1000 |
| token 总数 | 150439 | 133623 |
| 平均 token 长度 | 150.439 | 133.623 |
| 未知 token 数 | 0 | 0 |

在这份 tokenizer 和语料下，GraphBPE 的 SAFE 文本平均多 16.816 个 token，约长 12.58%。

## 9. 常见问题

### CMD 报 `unrecognized arguments: ` 后跟反引号

这是把 PowerShell 换行符粘贴到了 CMD。使用本文提供的单行命令即可；CMD 若要手动换行，应使用 `^`，不要使用反引号。

### 找不到 GraphBPE 词表

确认 `--vocab-prefix` 指向 `E:\GraphBpe+SAFE\SAFE-main\safe\data\tokenizer\vocab3000ring300\pretrain-token`，不带扩展名，并确认三个词表文件均存在。

### `cannot import name 'attr_as' from 'safe.utils'`

工具函数现已统一到 `safe\utils\__init__.py`，避免同名的 `safe\utils.py` 被 `safe\utils\` 包遮蔽。`attr_as`、`MolSlicer` 等原有 API 均通过 `safe.utils` 访问，GraphBPE 的 `safe.utils.token` 子模块也保留在同一包中。复制项目时请保留完整的 `safe\utils\` 文件夹。

### 切分进度结束但有失败记录

查看对应输出目录的 `failed.jsonl` 和 `summary.json`。非法或空 SMILES、带原子映射或通配原子的输入、严格模式下切分边界发生变化、分区或重建检查失败等，都会导致记录失败。不能仅凭进度到达 100% 判断所有分子都成功。

## 10. 本说明的核对范围

本文依据当前本地脚本与数据路径编写；两个批处理脚本的 `--help` 均已成功执行，输入表头、三个 GraphBPE 词表文件、自动保存比较结果的代码及已有结果文件均已核对。本次文档编写未重新运行批量切分，也未覆盖原有切分结果。

## 从 GitHub 下载完整数据（Git LFS）

本仓库使用 Git LFS 保存两份 ChEMBL 原始文件，以及 GraphBPE、BRICS 各自的 `records.jsonl` 和 `safe.txt`，共 6 个文件，约 1.70 GB。其余代码、词表、摘要和比较结果使用普通 Git 保存。

安装 Git LFS 后，在 CMD / Anaconda Prompt 中执行：

```cmd
git lfs install
git clone https://github.com/zeyuanyu195-design/molecular-tokenizer.git
cd molecular-tokenizer
git lfs pull
git lfs ls-files
```

如果已经克隆仓库，在仓库目录中运行 `git lfs pull` 即可补齐数据。建议使用 Git + Git LFS 获取项目；网页上的文件预览或 ZIP 下载是否包含完整 LFS 对象取决于仓库设置。

LFS 下载后的文件仍在 `safe/data/` 原目录中，批处理命令无需改变文件格式。本文前面的绝对路径针对原始本机目录；如果克隆到了其他位置，请将命令中的项目根目录替换为实际克隆位置。原始数据和切分结果可在本地继续使用，勿把只有几行内容的 LFS 指针文件当成完整数据。

重复上传修改后的大文件会增加 LFS 存储用量；下载完整数据会消耗仓库所有者的 LFS 下载额度。GitHub 的用量和计费说明见 https://docs.github.com/en/billing/concepts/product-billing/git-lfs 。
