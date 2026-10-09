# MolecularTokenizer

一个代码库，分别选择**分子切分方法**和**编码表示**。当前版本为 0.2.0。

| fragmentation | representation | 实现与边界 |
|---|---|---|
| `brics` | `safe` | BRICS 切分 + SAFE 字符串 + 序列 BPE |
| `npe` | `safe` | NPE 学习原子分组 + 保留原始属性的 SAFE 字符串 + 序列 BPE，主要新增功能 |
| `npe` | `demodiff` | DemoDiff 原有 NPE 节点/连接表示，保留严格保真检查 |
| `brics` | `demodiff` | 固定 BRICS 分组，收集片段词表，复用 DemoDiff 节点/连接编解码 |

SAFE 表示不需要在每条编码旁保存额外边表，但连接信息仍包含在字符串中，不能视为零开销。
BRICS 不训练切分规则；SAFE 组合训练序列词表，BRICS+DemoDiff 收集节点/边词表。
本项目不包含分子生成神经网络。

## 安装与测试

推荐 Python 3.11。Windows PowerShell：

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -e ".[test]"
.\.venv\Scripts\python.exe -m pytest -q
```

Linux/macOS 将 Python 路径换成 `.venv/bin/python`。安装固定提交的 SAFE 需要 Git 和网络。
核心验证环境为 RDKit 2025.9.6、tokenizers 0.23.2；依赖在 pyproject.toml 中固定。
旧 requirements.lock.txt 是本机历史环境清单，包含本地路径，不用于新环境安装，也不纳入 Git。
GitHub Actions 已配置 Windows/Linux Python 3.11；远程任务需要推送后才会运行。

## 快速使用：NPE + SAFE

```python
from molecular_tokenizer import MolecularTokenizer

def main():
    model = MolecularTokenizer(
        fragmentation="npe",
        representation="safe",
        bpe_scope="fragment",
        strict=True,
    )
    corpus = ["CCO", "CCN", "CCC", "CCCl", "CCOC(=O)C", "c1ccccc1"]
    report = model.train(
        corpus,
        vocab_size=256,       # 序列 BPE 的目标词表大小
        motif_vocab_size=16,  # NPE 图单元的目标词表大小，与上面分开
        ring_vocab_size=1,
        num_workers=1,
    )
    encoded = model.encode("[13CH3]CO")
    print(model.to_safe("[13CH3]CO"))
    print(encoded.token_ids, encoded.connections)  # SAFE 的 connections 为 ()
    print(model.decode(encoded))                  # [13CH3]CO
    model.save("my_npe_safe_model")              # 目标目录必须不存在
    restored = MolecularTokenizer.load("my_npe_safe_model")
    assert restored.encode("[13CH3]CO") == encoded

if __name__ == "__main__":
    main()  # Windows 上的多进程训练需要 main 保护
```

`examples/train_options.py` 是不依赖 ChEMBL 或外部仓库的四种组合示例。

## 分词范围与预算

`bpe_scope="fragment"`（新接口默认）让 BPE 在 SAFE 点分块内跨原子合并；点号保持为独立单元。
不解析 token 为独立分子：一个 BPE token 可能只是合法字符串的一部分。解码时先拼回完整串。
`bpe_scope="lexical"` 保留旧的原子/语法词素边界。DemoDiff 表示不接受 bpe_scope。

SAFE 的 vocab_size 是序列词表目标（至少 95，含完整基础字符表与 UNK）；实际大小可能小于目标。
NPE+SAFE 的 motif_vocab_size 控制图单元预算，ring_vocab_size 控制初始环系预算，二者独立于序列预算。
NPE+DemoDiff 的 vocab_size 直接控制图节点词表；BRICS+DemoDiff 将 vocab_size 作为上限，若观察到的
片段数超过上限会明确报错，不会静默丢弃片段。BRICS 图单元不是保留 dummy 标签的 BRICS 片段袋，
而是原子分组转换为 DemoDiff 模板，并用显式连接位置连接。

## 复用已有 NPE 词表

```python
model = MolecularTokenizer(fragmentation="npe", representation="safe")
report = model.train(
    smiles_list,
    vocab_size=3000,
    npe_model="path/to/saved_npe_model_directory",
)
```

npe_model 必须是通用类保存的、具有清单和校验值的 NPE 模型目录，也可以是之前保存的 NPE+SAFE。
直接的 `.node/.edge/.ring` 文件先用 `MolecularTokenizer.from_graphbpe(prefix).save(directory)` 导入。
导入时检查源模型，保存新模型时将 NPE 状态一并保存；加载新模型不依赖原路径。
SAFE 输出没有逐分子边表；如果复用旧模型，模型目录中的 graph.edge 可能仍保留其历史边词表，
这些边词表不用于 NPE+SAFE 编解码。

## 命令行

输入文件每行一条 SMILES；空行也会作为无效样本报告，不会悄悄改变源行编号。

```text
python -m molecular_tokenizer train --input train.smi --output model_npe_safe --fragmentation npe --representation safe --vocab-size 3000 --motif-vocab-size 350 --ring-vocab-size 100 --num-workers 4 --skip-invalid
python -m molecular_tokenizer encode --model model_npe_safe --smiles "[13CH3]CO"
python -m molecular_tokenizer evaluate --model model_npe_safe --input test.smi --output evaluation.json
```

用 `--max-molecules 10000` 只读取前 10,000 行（不是随机采样）。CLI 使用严格模式，evaluate 会逐条
记录错误，成功样本的平均长度与失败数同时输出。train 写入输入 SHA256、选取规则和拒绝源行。
不要用在训练数据上测得的结果声称泛化性能。

NPE+SAFE 的 num_workers 用于 NPE 图训练及 SAFE 转换；转换工作进程获得同一份冻结的 NPE 词表。
BRICS+SAFE 支持并行字符串转换；BRICS+DemoDiff 当前要求 num_workers=1。
NPE 默认内存模式随训练集增长。全量 NPE+SAFE 可选择 `--npe-storage sqlite --work-dir training_runs/npe_checkpoint`，
将图状态、全局候选频率和倒排索引保存在磁盘，分批训练；中断后加 `--resume` 使用同一目录续训。
断点会验证数据内容、顺序及配置，不接受不同输入；图训练可以续训，未完成的序列 BPE 阶段会重新运行。

## 保真策略与已知限制

NPE+SAFE 只使用 NPE 的原子分组，不调用 DemoDiff 的解码器。SAFE 从原始分子上切分，保留原子属性；
映射过程不重新编号原始分子。环内键不切断，包括 NPE 词表外的芳香环，避免产生无法解析的断环片段。
SAFE 的 ignore_stereo=False 还会保护不安全的立体切分。因此最终 SAFE 块可能不同于原始 NPE 分组。

strict=True 默认执行编码—解码后的规范化异构 SMILES 比较，不一致直接拒绝；不是自动修复。
MoleculeEncoding 只含 ID、连接和模型指纹，不存原始 SMILES。解码必须使用同一模型。
strict=False 仅关闭往返检查，不能解释为支持无损处理。

DemoDiff 表示仍继承上游的电荷、同位素和立体信息限制，BRICS 分组不会自动修好该解码器。
冻结图词表遇到未知单元或未知边会报错，不随机替换。平行片段连接当前不支持。
SAFE 暂不支持开放 dummy/polymer 输入；增强 CXSMILES 立体组等需要专门处理的输入可能被拒绝。

## 验证结果（2026-10-09）

本地 Windows Python 3.11 下 100 项测试通过；覆盖四种组合、片段内跨原子合并、特殊化学属性、
模型保存加载、复用词表、文件损坏检测、失败训练保留旧状态、拒绝源行和 CLI。已成功构建 wheel。
初始版本的 GitHub Actions Windows/Linux 测试已通过（[运行记录](https://github.com/zeyuanyu195-design/molecular-tokenizer/actions/runs/37889695528)）。

在已有 ChEMBL 随机样本中按规范化 SMILES 去重，seed=20261009 划分为 1,600 个训练分子和 400 个
测试分子，重新训练两种 SAFE 组合。二者序列词表目标均为 1,024，NPE 另用 350 个 motif / 100 个 ring：

| 组合 | 测试往返一致 | 测试平均 token 数 | 中位数 |
|---|---:|---:|---:|
| BRICS + SAFE + 片段内 BPE | 400/400 | 22.2875 | 20 |
| NPE + SAFE + 片段内 BPE | 400/400 | 25.27 | 21 |

这轮验证证明代码可运行且在该样本上保真，没有证明 NPE+SAFE 压缩更优；也不代表全量数据或下游表现。
复现实验脚本为 examples/validate_composable.py，需要历史 ChEMBL sample_selection.json；
其输入与生成模型属于本地数据，未放入 Git。历史 0.1 版全量/抽样模型没有覆盖或重新训练。

## 兼容性、文件与 GitHub

`MolecularTokenizer("safe")` 保持原 BRICS+SAFE 词法配置；`MolecularTokenizer("npe")` 保持原
NPE+DemoDiff。旧 schema=1 模型继续可加载。新代码将配置和 NPE 状态纳入指纹。
序列编码中的 backend 标签仍为 safe，DemoDiff 图编码仍为 npe，具体切分选择保存在模型配置里。

核心代码：tokenizer.py 管生命周期；fragmentation.py 管原子分组；pipelines.py 组合配置；
backends.py 与 _vendor/ 管兼容的原后端；types.py 管结构化输出；__main__.py 提供 CLI。

源码、测试、参考词表 fixtures 和 CI 配置统一维护在
[molecular-tokenizer](https://github.com/zeyuanyu195-design/molecular-tokenizer)。虚拟环境、ChEMBL 数据、
训练模型、分析结果和日志不纳入 Git。GitHub Actions 会在推送和拉取请求时运行测试；实际状态以仓库 Actions 页面为准。
第三方出处及许可见 THIRD_PARTY_NOTICES.txt；旧 README_zh.txt 保留历史实验说明。

## 全量实验与结果页面

全量 ChEMBL 实验的资源评估、可复现命令和分析口径见 [RESEARCH.md](RESEARCH.md)。
资源测试的抽样模型不等同于全量训练结果。磁盘版在同一批 5,000 条 ChEMBL 数据上训练至 3,000 个 NPE 单元后，
与原内存版的完整模型指纹相同。

2026-10-10 已完成全量训练和全部 2,897,819 条记录的评估。两种方法均准确还原
2,897,778 条；22 条输入无效，另有相同的 19 条立体化学相关不一致。
在共同准确还原的分子上，NPE+SAFE 平均 19.662 个 token，BRICS+SAFE 为 19.942，
NPE 少约 1.40%；但 p95 长度为 40 对 37，说明长序列尾部并未改善。
每种方法的序列词表均为 3,000；NPE 另有 3,000 个子图单元，模型容量并不相等。
这是训练集内描述性分析，不是独立测试泛化或下游性能结论。

[完整图表、词表和分析报告](https://zeyuanyu195-design.github.io/molecular-tokenizer/)
的静态源文件位于 `docs/index.html`。逐行数据与汇总的核验记录及失败诊断随报告提供。
