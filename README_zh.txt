通用分子分词器：SAFE 与 NPE（0.1 版历史说明）
2026-10-09：0.2 版已增加独立 fragmentation/representation 配置及 NPE+SAFE。
当前 API、安装、测试与限制请以同目录 README.md 为准；本文件保留历史实验记录。

已实现的入口：molecular_tokenizer.MolecularTokenizer
统一方法：train、encode、decode、encode_many、save、load。
附加入口：from_graphbpe，用于导入现有 GraphBPE 的 .node/.edge/.ring 词表。

设计依据
SAFE 是 Sequential Attachment-based Fragment Embedding；先把分子分成片段，
用 SMILES 环闭合标记表达片段连接，形成一个字符串表示。
NPE 是 DemoDiff 的 Node Pair Encoding；按邻接子图出现频率迭代学习 motif，
保留 motif 之间的键类型与连接原子位置，形成图表示。
因此，两者共享分词器的生命周期，但 encode 的结果需要保留不同的拓扑信息。

train 的区别
SAFE：官方 SAFEConverter 转换 SMILES，官方 SAFESplitter 进行化学语法分块，
随后用 Hugging Face tokenizers 的 Rust BPE 学习序列词表。
它支持迭代器训练，不跨分子边界合并；BRICS 等碎片规则本身不通过 BPE 学习。
完整 ASCII 语法字符表用于处理训练中未见的同位素数字等，避免静默 UNK 丢失。
这里使用 SAFE 官方转换和拆分，但不是 SAFE-GPT 的预训练 tokenizer/vocab。
NPE：复用 DemoDiff graphbpe.py 的环初始化与邻接 motif 频率合并，
之后遍历训练语料建立边词表；它仍把语料放入内存并重复访问分子。
SAFE 的 vocab_size 是包含基础语法字符和 UNK 的词表目标；
NPE 的 vocab_size 是节点词表目标，ring_vocab_size 控制初始环词表。
两者都可能在耗尽可学习候选时返回小于目标值的实际词表。

encode 的区别
SAFE：SMILES -> SAFE 字符串 -> 序列 token IDs。
NPE：SMILES -> motif 节点 IDs + 显式 Connection 列表。
Connection 保存 source、target、bond_type、source_attachment、target_attachment。
所有连接只保存一次，source < target；键类型约定 single=0、double=1、triple=2。
NPE 的边词表默认冻结；即使 strict=False，也拒绝未知连接，不会丢弃未知键。

decode 的区别
SAFE：token IDs -> 完整 SAFE 字符串 -> 官方转换器 -> 规范化 SMILES。
NPE：节点 IDs + 键类型与连接位置 -> 原 GraphBPE decode -> 规范化 SMILES。
NPE 的 node IDs 单独不足以解码一个分子。
不会从原始 SMILES、原始分子对象或缓存重建，输出载荷中没有这些信息。
SAFE 解码不自动修复不完整的生成字符串，不自动删除 dummy 原子。
NPE 不随机替换未知 motif。

公共载荷
MoleculeEncoding.backend：safe 或 npe。
MoleculeEncoding.tokenizer_id：词表和后端配置的 SHA-256 身份。
MoleculeEncoding.token_ids：SAFE 序列 ID，或 NPE 节点 ID。
MoleculeEncoding.connections：NPE 连接；SAFE 时为空。
MoleculeEncoding.fidelity_checked：是否在 encode 时完成严格分子往返核对。
to_dict/from_dict 可用于 JSON 序列化，不要求 PyTorch 或模型权重。
同一公共类不意味着序列模型和图模型可以直接接收对方的输入。
后续 Transformer/DiT 的模型输入整理仍需按 backend 处理。

最小使用示例
from molecular_tokenizer import MolecularTokenizer

corpus = ['CCO', 'CCN', 'CCCl', 'c1ccccc1', 'Cc1ccccc1']
safe = MolecularTokenizer('safe')
safe.train(iter(corpus), vocab_size=256)
encoded = safe.encode('CCCl')
smiles = safe.decode(encoded)

npe = MolecularTokenizer('npe')
npe.train(corpus, vocab_size=12, ring_vocab_size=1, num_workers=2)
encoded_graph = npe.encode('CCCl')
smiles = npe.decode(encoded_graph)

在 Windows 中，包含 NPE train 的入口脚本必须写成：
if __name__ == '__main__':
    main()
完整可运行例子在 examples/demo.py。

导入已经部署的 NPE 词表
from pathlib import Path
prefix = Path(r'E:\codex\2026-10-07\prepare-all-the-smiles-strings-from\outputs\DemoDiff\data\tokenizer\vocab3000ring300\pretrain-token')
npe = MolecularTokenizer.from_graphbpe(prefix)

模型保存和加载
model_dir = safe.save('my_new_safe_tokenizer')
loaded = MolecularTokenizer.load(model_dir)
assert loaded.decode(encoded) == safe.decode(encoded)
save 需要一个尚不存在的目录；它先完成写入再原子移动目录，不覆盖现有模型。
tokenizer.json 保存后端、配置、版本与各个必要词表文件的 SHA-256。
load 会核对文件校验和和词表身份，拒绝损坏或混用的模型。
重新训练先使用临时后端对象，失败时原模型仍可用。
单实例公开方法使用锁串行保护；不要直接修改 backend 内部词表。

直接运行（PowerShell，从项目目录）
项目目录：E:\codex\2026-10-07\prepare-all-the-smiles-strings-from\outputs\generic_tokenizer
& .\.venv\Scripts\python.exe -X utf8 -m pytest -q
& .\run_demo.ps1 --limit 100
演示读取之前准备的 ChEMBL .smi 文件前 100 条，训练 SAFE 小词表，
训练 NPE 玩具词表并导入 DemoDiff 预训练词表，再生成序列化、重载与解码报告。
每次演示创建独立 validation/demo-* 文件夹；summary 保存在
validation/latest_demo_summary.json。

安装与版本
本项目已单独部署在 .venv 中，原 minbpe 和 DemoDiff 环境不受影响。
运行版本：Python 3.11.7，RDKit 2025.9.6，tokenizers 0.23.2，NumPy 2.4.6。
SAFE 使用固定源码提交 d162d23e3e8c8a77b54c559e88a35d08b3b1d1b5。
初次在其他位置安装：Python 3.11 创建 .venv，随后 pip install -e '.[test]'。
若需复现已安装依赖，先 pip install -r requirements.lock.txt，
再 pip install --no-deps -e .。
不同 RDKit 版本可能改变规范化、片段编号或图身份判断，应重新验证模型。

验证结果
34 项测试通过，包含两个后端的往返、序列化、模型保存加载、
词表不变性、不同词表拒绝、训练失败回滚、损坏模型拒绝、
无环语料训练、空边/环词表重载、Windows 多进程与特殊结构。
ChEMBL 37 全量 SAFE 训练接受 2,897,797/2,897,819 条，记录 22 条拒绝项，
目标词表 3000，实际序列词表 1176。
NPE 全量构图在本机内存限制下不可行；最终使用固定 seed=42 的 200,000 条
均匀蓄水池样本训练，目标节点词表 350（其中初始 ring_vocab_size=300）。
训练清单会记录完整来源数据哈希、实际样本行和采样配置；不将此结果称为全量 NPE 词表。

复现当前 ChEMBL 词表
在 outputs/generic_tokenizer 目录运行：

python examples/train_chembl37.py safe --vocab-size 3000 --num-workers 12
python examples/train_chembl37.py npe --vocab-size 350 --num-workers 8 --max-rows 200000 --seed 42
python examples/compare_vocabularies.py

SAFE 词表在 trained_models/chembl37_safe_vocab3000/vocabulary.tsv；
NPE 词表在 trained_models/chembl37_npe_vocab350ring300/vocabulary.tsv。
NPE 的 training_selection.tsv 保存所选样本的源行号和 ChEMBL ID。

严格性与当前范围
strict=True 是默认设置；encode 比较 RDKit 规范化的含立体信息分子身份，
不要求重建字符串与输入文字逐字相同。盐、同位素和立体信息参与身份判断。
现有 NPE 对部分原子电荷、同位素及跨 motif 的立体信息仍可能丢失；
本实现检测并拒绝这些情况，没有通过原始 SMILES 缓存绕过解码。
strict=False 可明确关闭往返核对，此时 fidelity_checked=False；
这不能保证结构信息保留。它仍拒绝未知 motif、未知连接和非法图载荷。
现有 GraphBPE decoder 不支持两 motif 间的并行多键，适配层显式拒绝。
SAFE 的开放 dummy/polymer 连接点需要另一种 SAFE fragment 工作流；
当前 SMILES 完整分子接口对这些输入显式拒绝。
encode_many 为惰性迭代器，任何失败都抛出异常；做数据清洗时按条捕获异常并
记录拒绝原因，不要默认忽略失败样本。

上游来源
SAFE：https://github.com/datamol-io/safe/tree/d162d23e3e8c8a77b54c559e88a35d08b3b1d1b5
SAFE 论文：https://doi.org/10.1039/D4DD00019F
NPE/GraphBPE：https://github.com/liugangcode/DemoDiff/blob/0acb85ef446c175a3e3a60bbad43de97d97ca76f/downstream/graphbpe.py
NPE 论文：https://arxiv.org/html/2510.08744v1#S3.SS1
GraphBPE 的 MIT 许可证与五项局部兼容性修改记录保存在 _vendor 目录。
原先部署的两个仓库仍保持未修改状态。
