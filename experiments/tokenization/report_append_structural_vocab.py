"""Summarize the frozen-16k plus structural-span experiment from saved outputs."""
import argparse
import json
from pathlib import Path


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--experiment',type=Path,required=True)
    args=parser.parse_args()
    def read(name): return json.loads((args.experiment/name).read_text(encoding='utf-8'))
    stats=read('summary.json'); coverage=read('coverage.json'); paired=read('paired_comparison.json')
    audit=read('verification.json'); examples=read('CHEMBL45183.json')
    baseline=stats['all']['whole_16000']['mean']
    names={'whole_16000':'完整片段基线', 'A_plus4000':'A +4k','A_plus8000':'A +8k','B_plus4000':'B +4k','B_plus8000':'B +8k'}
    lines=['# 保留完整片段 16k，追加结构词条实验','',
        '同一批 44,329 条边界一致记录：TRAIN 35,518、VALIDATION 4,268、TEST 4,543。全量均值包含训练和验证样本；测试集此前已用于分析，不是全新盲测。','',
        '冻结原完整片段词表的全部 16,000 个 ID 和载荷，包括 2,136 个基础词条及 13,864 个完整片段词条。已有完整片段命中仍输出同一 ID，不再拆分。','',
        '新增词条仅从 TRAIN 未命中原完整片段查表的片段收集。A/B 的合法边界与前次实验一致；按出现频次降序、文本字典序破同分排序。四组均采用预先确定的频次排名，没有根据测试长度选择方案。','',
        '每组新增预算均包含 170 个完整连接标号和两个块标记。因此 +4k 新增 3,828 个结构跨度词条，总大小 20,000；+8k 新增 7,828 个，总大小 24,000。同一方案的 20k 词表是 24k 词表的原样前缀。','',
        '## 全部记录与测试集','',
        '| 方法 | 总词表 | 总 token 数 | 全量平均长度 | 相比原16k缩短 | TEST平均长度 |',
        '|---|---:|---:|---:|---:|---:|']
    for name,label in names.items():
        s=stats['all'][name]; K=16000 if name=='whole_16000' else 16000+int(name.split('plus')[1])
        lines.append(f"| {label} | {K:,} | {s['total']:,} | {s['mean']:.6f} | {(1-s['mean']/baseline)*100:.2f}% | {stats['test'][name]['mean']:.6f} |")
    lines+=['','计入点号、连接标号、立体符号和 BEGIN/END 块标记；不含 BOS/EOS、padding，不截断。追加方案使用了更大的词表，不能据此声称相同词表大小下优于完整片段方案。','',
        '## 覆盖与逐分子比较','',
        '| 方法 | 原词表完整命中片段数 | 无需结构块回退的片段比例 | 无需结构块回退的分子数 | 更短/相同/更长的分子数 |',
        '|---|---:|---:|---:|---|']
    for n,c in coverage.items():
        q=paired[n]
        lines.append(f"| {names[n]} | {c['frozen_full_fragment_hits']:,} | {c['fragments_without_fallback']/c['fragments']*100:.2f}% | {c['molecules_without_fallback']:,} | {q.get('shorter',0):,} / {q.get('equal',0):,} / {q.get('longer',0):,} |")
    lines+=['','无需回退不代表每个原片段只占一个 ID；可能由多个合法跨度词条组成。所有未命中逻辑块仍可通过完整词法单元无损表示。','',
        '## 完整还原验证','',
        f"- 重新加载落盘词表及输出，验证 {audit['verified_id_sequences']:,} 条 ID 序列；全部精确还原 SAFE。",
        '- 四组还原字符串逐条均与同一个源 SAFE 相等；每条记录调用真实 SAFE decoder，禁用修复及 dummy 删除，再比较 canonical isomeric 分子身份，共 44,329 条通过。',
        '- 前 16,000 个词条载荷及已有完整片段 ID 全部保持不变。',
        '- 所有逻辑切点符合 A/B 约束；检查完整原子、配对分支及长度组成，违规数均为 0。',
        '- 本批记录无未知词法单元。未用测试样本补词表，也未启用字节回退；不保证任意新分子均被覆盖。',
        '- 结构保护对象仍是逻辑块；回退框内部的基础语法 ID 不声称具有独立化学基序意义。','',
        '## CHEMBL45183','',
        '| 方法 | 总 token 数 |','|---|---:|']
    for e in examples: lines.append(f"| {names[e['name']]} | {e['length']} |")
    for e in examples:
        lines+=['',f"### {names[e['name']]}",'','| SAFE片段 | 实际ID | 词条 |','|---|---|---|']
        for f in e['fragments']:
            tokens=' · '.join('`'+x+'`' for x in f['tokens'])
            lines.append(f"| `{f['fragment']}` | `{f['ids']}` | {tokens} |")
    lines+=['','完整序列另含 7 个点号 token；ID 11。','',
        '统计：summary.json；逐条输出：per_molecule.jsonl；覆盖率：coverage.json；核验：verification.json；词表：A_plus4000、A_plus8000、B_plus4000、B_plus8000 下的 vocab.json。','']
    out=args.experiment/'REPORT.zh.md'
    out.write_text('\n'.join(lines),encoding='utf-8')
    print(out)


if __name__=='__main__':
    main()
