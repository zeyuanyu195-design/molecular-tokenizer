# 完整语法块与受保护子结构词表实验

## 结论口径

同一批 44,329 条 GraphBPE 与 SAFE 边界一致记录；训练 35,518，验证 4,268，测试 4,543。仅 TRAIN 学词表，用验证集在频次/估计节省量两种排名间选择，再评估测试集。完整 SAFE 文本和原始点号分片均不改变。

A 和 B 均保留完整片段优先匹配；未命中时以动态规划寻找词表下最短序列。词条匹配只能在合法逻辑块的边界间进行；匹配失败的块使用 BEGIN/END 包围的完整原子、完整连接标号及语法单元表示。

**保护单位是逻辑子结构，不是每一个模型 token。** 回退框内的括号、键和标号是结构语法 token。当前方案不宣称每个模型 token 都是独立完整基序，也不宣称所有连通子图均具有特定生物化学功能。

## A 与 B 的明确规则

- A：方括号原子、多字符原子、`%12`/`%(100)` 等连接标号均不可拆分。一个顶层原子连同它的完整分支和环标号组成基本块；进入下一原子的键符号随下一块保存。候选由连续基本块组成，因此不会出现 `CCC(C` 或 `)C45` 这样的内容词条。
- B：在 A 基础上，从完整 SAFE 解析分子图，对原 SAFE 片段内部的环原子组、`[CX3](=[OX1])[NX3]` 酰胺组、`[CX3](=[OX1])[OX2]` 羧酸/酯组禁止内部切分。重叠保护区间合并。检查每块诱导子图连通且相邻块有图边，因而每个连续候选跨度也连通。
- 已跨原始点号片段的环或官能团只记录数量，不重新合并原片段；不能将 B 理解为保护整个分子中所有官能团。
- 候选最多 16 个基本块、32 个原子、128 个完整词法单元；这只限制新增子串候选，不限制完整片段词条或大块的无损兜底。

## 词表预算和对照限制

- 原基础层 2,136 保留。补入 170 个训练集完整连接标号词条和两个块标记，总基础层 M=2,308。原基础层已有的完整稀有原子（如 `[Mg+2]`）继续使用，没有从测试集添加词条。
- 16k：完整片段仍为 6,932，新增结构跨度词条 6,760。32k：完整片段仍为 7,466，结构跨度词条 22,226。总大小严格为 16,000 / 32,000。
- 原基础层含不完整组合词条，但编码只在合法结构块边界匹配，不能输出跨半个原子/分支的内容词条；兜底框内只使用完整词法单元。
- A 与 B 的算法和预算一致，主要区别是结构保护产生的合法切分边界。

## 全部 44,329 条的平均 token 长度

| 方法 | 16k | 32k |
|---|---:|---:|
| A：完整语法块 | 26.936588 | 24.891065 |
| B：保护环系及选定官能团 | 35.162016 | 28.663358 |

计入连接编号、点号、键、立体符号、兜底块标记；不含 BOS/EOS/padding，不截断。全量结果包含训练和验证样本，不代替测试集泛化结果。

## 测试集 4,543 条

| 方法 | 16k | 32k |
|---|---:|---:|
| A | 28.464891 | 26.302884 |
| B | 37.098833 | 31.513978 |

## 验证集选择

| 方法 | K | 选中的排名方式 | 验证均值 |
|---|---:|---|---:|
| A | 16000 | frequency | 26.582474 |
| A | 32000 | frequency | 24.575914 |
| B | 16000 | frequency | 35.479381 |
| B | 32000 | frequency | 29.909091 |

频次排名按训练跨度出现次数降序；gain 排名按训练出现次数乘以相对基础块表达的估计节省量。跨度可重叠，gain 不是实际全语料边际收益。两个候选均在验证集比较，选择结果先落盘再评估测试集。

## 完整覆盖和结构检查

| 方法 | ID/SAFE/分子还原 | 使用兜底的分子数 | 兜底块数 | 标记 token 数 |
|---|---:|---:|---:|---:|
| A_K16000_frequency | 44,329/44,329 | 18,341 | 54,101 | 108,202 |
| A_K32000_frequency | 44,329/44,329 | 14,560 | 44,106 | 88,212 |
| B_K16000_frequency | 44,329/44,329 | 31,455 | 73,782 | 147,564 |
| B_K32000_frequency | 44,329/44,329 | 22,220 | 46,425 | 92,850 |

所有方法的 ID 序列精确还原相同原 SAFE；每条记录再经 SAFE decoder（禁用修复、禁用 dummy 删除）与原分子的 canonical isomeric SMILES 比较。独立重新读取落盘 IDs，验证合法逻辑切点、完整原子和括号、长度求和及源文件哈希。

本批数据未出现原基础词表加 TRAIN 词法单元仍不覆盖的输入。没有启用破坏原子单元的字节兜底；对于其他新数据中完全未见的原子写法/编号，当前实现会显式失败，不能宣称任意分子全覆盖。

保护组计数（按记录出现次数，可能重叠）：

```json
{
  "protected_ring": 110433,
  "A_blocks": 1137491,
  "B_blocks": 616661,
  "protected_amide": 35244,
  "cross_fragment_amide": 9520,
  "cross_fragment_ring": 34960,
  "protected_carboxylic_acid_or_ester": 15764,
  "cross_fragment_carboxylic_acid_or_ester": 2486
}
```

## CHEMBL45183 的真实输出

### A_K16000_frequency：21 个 token

| 原 SAFE 片段 | 允许的基本逻辑块 | 最终模型 token 内容 |
|---|---|---|
| `C123CCC(C1)C45C26` | `C123` · `C` · `C` · `C(C1)` · `C45` · `C26` | `C123` · `CCC(C1)` · `C45` · `C26` |
| `C7(=O)C8NC(=O)[C@H]9N` | `C7(=O)` · `C8` · `N` · `C(=O)` · `[C@H]9` · `N` | `C7(=O)` · `C8N` · `C(=O)[C@H]9N` |
| `C8c1cccs1` | `C8` · `c1` · `c` · `c` · `c` · `s1` | `C8c1ccc` · `s1` |
| `C9C(=O)O` | `C9` · `C(=O)` · `O` | `C9C(=O)O` |
| `O67` | `O67` | `O67` |
| `C3` | `C3` | `C3` |
| `C4` | `C4` | `C4` |
| `C5` | `C5` | `C5` |

### A_K32000_frequency：21 个 token

| 原 SAFE 片段 | 允许的基本逻辑块 | 最终模型 token 内容 |
|---|---|---|
| `C123CCC(C1)C45C26` | `C123` · `C` · `C` · `C(C1)` · `C45` · `C26` | `C123CC` · `C(C1)` · `C45` · `C26` |
| `C7(=O)C8NC(=O)[C@H]9N` | `C7(=O)` · `C8` · `N` · `C(=O)` · `[C@H]9` · `N` | `C7(=O)` · `C8NC(=O)` · `[C@H]9N` |
| `C8c1cccs1` | `C8` · `c1` · `c` · `c` · `c` · `s1` | `C8c1ccc` · `s1` |
| `C9C(=O)O` | `C9` · `C(=O)` · `O` | `C9C(=O)O` |
| `O67` | `O67` | `O67` |
| `C3` | `C3` | `C3` |
| `C4` | `C4` | `C4` |
| `C5` | `C5` | `C5` |

### B_K16000_frequency：35 个 token

| 原 SAFE 片段 | 允许的基本逻辑块 | 最终模型 token 内容 |
|---|---|---|
| `C123CCC(C1)C45C26` | `C123CCC(C1)C45C26` | `<BLOCK_BEGIN>` · `C` · `1` · `2` · `3` · `C` · `C` · `C` · `(` · `C` · `1` · `)` · `C` · `4` · `5` · `C` · `2` · `6` · `<BLOCK_END>` |
| `C7(=O)C8NC(=O)[C@H]9N` | `C7(=O)` · `C8` · `NC(=O)` · `[C@H]9` · `N` | `C7(=O)C8` · `NC(=O)[C@H]9N` |
| `C8c1cccs1` | `C8` · `c1cccs1` | `C8` · `c1cccs1` |
| `C9C(=O)O` | `C9` · `C(=O)O` | `C9C(=O)O` |
| `O67` | `O67` | `O67` |
| `C3` | `C3` | `C3` |
| `C4` | `C4` | `C4` |
| `C5` | `C5` | `C5` |

### B_K32000_frequency：17 个 token

| 原 SAFE 片段 | 允许的基本逻辑块 | 最终模型 token 内容 |
|---|---|---|
| `C123CCC(C1)C45C26` | `C123CCC(C1)C45C26` | `C123CCC(C1)C45C26` |
| `C7(=O)C8NC(=O)[C@H]9N` | `C7(=O)` · `C8` · `NC(=O)` · `[C@H]9` · `N` | `C7(=O)C8` · `NC(=O)[C@H]9N` |
| `C8c1cccs1` | `C8` · `c1cccs1` | `C8` · `c1cccs1` |
| `C9C(=O)O` | `C9` · `C(=O)O` | `C9C(=O)O` |
| `O67` | `O67` | `O67` |
| `C3` | `C3` | `C3` |
| `C4` | `C4` | `C4` |
| `C5` | `C5` | `C5` |

## 核心伪代码

```text
TRAIN(mode, K):
    preserve the previous whole-fragment ranking
    form complete lexical/branch blocks from TRAIN SAFE strings
    if mode == B: coalesce blocks across in-fragment protected groups
    enumerate valid contiguous spans; count TRAIN occurrences
    build frequency and estimated-gain vocabularies within K
    select ranking using validation total token length; freeze vocabulary

ENCODE(safe, mode, vocab):
    parse full SAFE to establish atom/block membership
    for original SAFE fragment:
        if whole fragment is in vocabulary: emit its ID
        else:
            construct legal blocks for mode
            use DP to minimize token count over dictionary span matches
            fallback edge: BEGIN + complete atomic lexeme IDs + END
    preserve dots between original fragments

DECODE(ids, vocab):
    expand content IDs and framed primitive IDs
    remove framing markers to obtain the exact original SAFE
    decode complete SAFE and validate molecular stereochemistry
```

## 文件

- 每个选定方案目录的 vocab.json：完整 ID 词表和解码信息。
- per_molecule.jsonl：44,329 行，含四组实际 ID、长度、兜底次数和原基线长度。
- partition_cache.jsonl：每条记录的 A/B 合法逻辑块和保护组统计。
- CHEMBL45183.json：上述示例的实际 IDs 和内容。
- summary.json / verification.json：各数据划分结果和独立检查。
