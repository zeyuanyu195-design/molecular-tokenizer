"""Independent persisted-ID / logical-cut audit and Chinese experiment report."""
import argparse
from collections import Counter
from functools import lru_cache
import json
from pathlib import Path

from evaluate_structural_span_vocab import SpanCodec,lex,ATOM,syntactic_units
from evaluate_bounded_safe_vocab import dump,sha


@lru_cache(maxsize=100000)
def check_content(text):
    tokens=lex(text)
    if not any(ATOM.fullmatch(t) for t in tokens):
        raise ValueError(f'Content token has no atom: {text}')
    depth=0
    for t in tokens:
        depth+=(t=='(')-(t==')')
        if depth<0:
            raise ValueError(f'Content token begins inside a branch: {text}')
    if depth:
        raise ValueError(f'Content token has unmatched branch: {text}')
    syntactic_units(text)
    return True


def logical_outputs(codec,ids):
    pieces=[]
    i=0
    while i<len(ids):
        idx=ids[i]
        if idx==codec.dot:
            pieces.append(('.',False))
        elif idx==codec.begin:
            start=i
            i+=1
            while i<len(ids) and ids[i]!=codec.end:
                i+=1
            if i==len(ids):
                raise ValueError('Unclosed frame in persisted IDs')
            pieces.append((codec.decode(ids[start:i+1]),True))
        else:
            pieces.append((codec.payloads[idx].decode(),False))
        i+=1
    return pieces


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--experiment',type=Path,required=True)
    p.add_argument('--records',type=Path,required=True)
    p.add_argument('--previous',type=Path,required=True)
    args=p.parse_args()
    folder=args.experiment
    selected=json.loads((folder/'selection_before_test.json').read_text())['selected']
    codecs={name:SpanCodec(json.loads((folder/name/'vocab.json').read_text())) for name in selected.values()}
    records={r['record_index']:r for r in map(json.loads,args.records.read_text().splitlines())}
    partitions={r['record_index']:r['partitions'] for r in map(json.loads,(folder/'partition_cache.jsonl').read_text().splitlines())}
    totals=Counter()
    frames=Counter()
    framed_molecules=Counter()
    literal_units=Counter()
    primitive_cost=Counter()
    dot_cost=Counter()
    seen=set()
    for line in (folder/'per_molecule.jsonl').open():
        r=json.loads(line)
        idx=r['record_index']
        assert idx not in seen
        seen.add(idx)
        assert r['cohort']==records[idx]['cohort']
        for name,value in r['methods'].items():
            c=codecs[name]
            assert len(c.payloads)==c.config['K']
            assert c.decode(value['ids'])==records[idx]['safe']
            assert value['length']==len(value['ids'])
            assert value['text_roundtrip'] and value['molecular_roundtrip']
            parts=partitions[idx][c.config['mode']]
            boundaries=[]
            for units in parts:
                end=0
                allowed={0}
                for u in units:
                    end+=len(u)
                    allowed.add(end)
                boundaries.append(allowed)
            fi=0
            pos=0
            fcount=0
            for text,is_frame in logical_outputs(c,value['ids']):
                if text=='.':
                    assert pos==max(boundaries[fi])
                    fi+=1
                    pos=0
                    dot_cost[name]+=1
                    continue
                check_content(text)
                pos+=len(text)
                # Every logical cut must be at an A/B-approved structural cut.
                # In B, these cuts cannot split any in-fragment protected group.
                assert pos in boundaries[fi]
                if is_frame:
                    fcount+=1
                    primitive_cost[name]+=len(lex(text))
                else:
                    literal_units[name]+=1
            assert fi==len(parts)-1 and pos==max(boundaries[fi])
            assert fcount==value['fallback_blocks']
            frames[name]+=fcount
            framed_molecules[name]+=fcount>0
            totals[name]+=value['length']
    assert seen==set(records)==set(partitions) and len(seen)==44329
    summary=json.loads((folder/'summary.json').read_text())
    for name,c in codecs.items():
        assert summary['all'][name]['total']==totals[name]
        assert totals[name]==literal_units[name]+primitive_cost[name]+2*frames[name]+dot_cost[name]
        # New learned spans / whole entries all have complete lexical structure.
        for payload in c.payloads[c.config['M']:]:
            check_content(payload.decode())
    provenance=json.loads((folder/'provenance.json').read_text())
    for f,expected in provenance['sources'].items():
        assert sha(f)==expected,f
    audit=dict(records=len(seen),verified_sequences=len(seen)*len(codecs),
        exact_id_safe_roundtrips=True,logical_cut_violations=0,content_syntax_violations=0,
        source_hashes_verified=True,per_molecule_sha256=sha(folder/'per_molecule.jsonl'),
        length_components={name:dict(literal_units=literal_units[name],primitive_ids=primitive_cost[name],
            frame_markers=2*frames[name],dots=dot_cost[name],fallback_blocks=frames[name],
            molecules_using_fallback=framed_molecules[name]) for name in codecs})
    dump(folder/'verification.json',audit)
    old=json.loads((args.previous/'summary.json').read_text())
    stats=json.loads((folder/'partition_audit.json').read_text())
    examples=json.loads((folder/'CHEMBL45183.json').read_text())
    grid={r['name']:r for r in json.loads((folder/'validation_grid.json').read_text())}
    md=['# 完整语法块与受保护子结构词表实验','',
        '## 结论口径','',
        '同一批 44,329 条 GraphBPE 与 SAFE 边界一致记录；训练 35,518，验证 4,268，测试 4,543。仅 TRAIN 学词表，用验证集在频次/估计节省量两种排名间选择，再评估测试集。完整 SAFE 文本和原始点号分片均不改变。',
        '',
        'A 和 B 均保留完整片段优先匹配；未命中时以动态规划寻找词表下最短序列。词条匹配只能在合法逻辑块的边界间进行；匹配失败的块使用 BEGIN/END 包围的完整原子、完整连接标号及语法单元表示。',
        '',
        '**保护单位是逻辑子结构，不是每一个模型 token。** 回退框内的括号、键和标号是结构语法 token。当前方案不宣称每个模型 token 都是独立完整基序，也不宣称所有连通子图均具有特定生物化学功能。',
        '',
        '## A 与 B 的明确规则','',
        '- A：方括号原子、多字符原子、`%12`/`%(100)` 等连接标号均不可拆分。一个顶层原子连同它的完整分支和环标号组成基本块；进入下一原子的键符号随下一块保存。候选由连续基本块组成，因此不会出现 `CCC(C` 或 `)C45` 这样的内容词条。',
        '- B：在 A 基础上，从完整 SAFE 解析分子图，对原 SAFE 片段内部的环原子组、`[CX3](=[OX1])[NX3]` 酰胺组、`[CX3](=[OX1])[OX2]` 羧酸/酯组禁止内部切分。重叠保护区间合并。检查每块诱导子图连通且相邻块有图边，因而每个连续候选跨度也连通。',
        '- 已跨原始点号片段的环或官能团只记录数量，不重新合并原片段；不能将 B 理解为保护整个分子中所有官能团。',
        '- 候选最多 16 个基本块、32 个原子、128 个完整词法单元；这只限制新增子串候选，不限制完整片段词条或大块的无损兜底。',
        '',
        '## 词表预算和对照限制','',
        '- 原基础层 2,136 保留。补入 170 个训练集完整连接标号词条和两个块标记，总基础层 M=2,308。原基础层已有的完整稀有原子（如 `[Mg+2]`）继续使用，没有从测试集添加词条。',
        '- 16k：完整片段仍为 6,932，新增结构跨度词条 6,760。32k：完整片段仍为 7,466，结构跨度词条 22,226。总大小严格为 16,000 / 32,000。',
        '- 原基础层含不完整组合词条，但编码只在合法结构块边界匹配，不能输出跨半个原子/分支的内容词条；兜底框内只使用完整词法单元。',
        '- 与旧 BPE 的比较同时改变了候选构建、匹配算法、基础单元和兜底标记，不能把全部长度差归因于一个括号约束。A 与 B 的算法和预算一致，主要区别是结构保护产生的合法切分边界。',
        '',
        '## 全部 44,329 条的平均 token 长度','',
        '| 方法 | 16k | 32k |','|---|---:|---:|',
        f"| 原混合 BPE | {old['all']['historical_16000']['mean']:.6f} | {old['all']['historical_32000']['mean']:.6f} |",
        f"| 仅禁止新增纯数字/符号组合 | {old['all']['constrained_K16000_whole50']['mean']:.6f} | {old['all']['constrained_K32000_whole25']['mean']:.6f} |"]
    for mode,label in [('A','A：完整语法块'),('B','B：保护环系及选定官能团')]:
        md.append(f"| {label} | {summary['all'][selected[mode+'_16000']]['mean']:.6f} | {summary['all'][selected[mode+'_32000']]['mean']:.6f} |")
    md+=['','计入连接编号、点号、键、立体符号、兜底块标记；不含 BOS/EOS/padding，不截断。全量结果包含训练和验证样本，不代替测试集泛化结果。','',
         '## 测试集 4,543 条','',
         '| 方法 | 16k | 32k |','|---|---:|---:|']
    for mode in ('A','B'):
        md.append(f"| {mode} | {summary['test'][selected[mode+'_16000']]['mean']:.6f} | {summary['test'][selected[mode+'_32000']]['mean']:.6f} |")
    md+=['','## 验证集选择','',
         '| 方法 | K | 选中的排名方式 | 验证均值 |','|---|---:|---|---:|']
    for key,name in selected.items():
        g=grid[name]
        md.append(f"| {g['mode']} | {g['K']} | {g['ranking']} | {g['validation']['mean']:.6f} |")
    md+=['','频次排名按训练跨度出现次数降序；gain 排名按训练出现次数乘以相对基础块表达的估计节省量。跨度可重叠，gain 不是实际全语料边际收益。两个候选均在验证集比较，选择结果先落盘再评估测试集。','',
         '## 完整覆盖和结构检查','',
         '| 方法 | ID/SAFE/分子还原 | 使用兜底的分子数 | 兜底块数 | 标记 token 数 |','|---|---:|---:|---:|---:|']
    for name in codecs:
        a=audit['length_components'][name]
        md.append(f"| {name} | 44,329/44,329 | {a['molecules_using_fallback']:,} | {a['fallback_blocks']:,} | {a['frame_markers']:,} |")
    md+=['','所有方法的 ID 序列精确还原相同原 SAFE；每条记录再经 SAFE decoder（禁用修复、禁用 dummy 删除）与原分子的 canonical isomeric SMILES 比较。独立重新读取落盘 IDs，验证合法逻辑切点、完整原子和括号、长度求和及源文件哈希。','',
         '本批数据未出现原基础词表加 TRAIN 词法单元仍不覆盖的输入。没有启用破坏原子单元的字节兜底；对于其他新数据中完全未见的原子写法/编号，当前实现会显式失败，不能宣称任意分子全覆盖。','',
         '保护组计数（按记录出现次数，可能重叠）：','',
         '```json',json.dumps(stats['counts'],ensure_ascii=False,indent=2),'```','',
         '## CHEMBL45183 的真实输出','']
    for example in examples:
        md += [f"### {example['name']}：{example['length']} 个 token",'',
               '| 原 SAFE 片段 | 允许的基本逻辑块 | 最终模型 token 内容 |','|---|---|---|']
        for f in example['fragments']:
            units=' · '.join('`'+t+'`' for t in f['logical_units'])
            tokens=' · '.join('`'+t+'`' for t in f['tokens'])
            md.append(f"| `{f['fragment']}` | {units} | {tokens} |")
        md.append('')
    md+=['## 核心伪代码','',
         '```text','TRAIN(mode, K):','    preserve the previous whole-fragment ranking',
         '    form complete lexical/branch blocks from TRAIN SAFE strings',
         '    if mode == B: coalesce blocks across in-fragment protected groups',
         '    enumerate valid contiguous spans; count TRAIN occurrences',
         '    build frequency and estimated-gain vocabularies within K',
         '    select ranking using validation total token length; freeze vocabulary','',
         'ENCODE(safe, mode, vocab):','    parse full SAFE to establish atom/block membership',
         '    for original SAFE fragment:','        if whole fragment is in vocabulary: emit its ID',
         '        else:','            construct legal blocks for mode',
         '            use DP to minimize token count over dictionary span matches',
         '            fallback edge: BEGIN + complete atomic lexeme IDs + END',
         '    preserve dots between original fragments','',
         'DECODE(ids, vocab):','    expand content IDs and framed primitive IDs',
         '    remove framing markers to obtain the exact original SAFE',
         '    decode complete SAFE and validate molecular stereochemistry','```','',
         '## 文件','',
         '- 每个选定方案目录的 vocab.json：完整 ID 词表和解码信息。',
         '- per_molecule.jsonl：44,329 行，含四组实际 ID、长度、兜底次数和原基线长度。',
         '- partition_cache.jsonl：每条记录的 A/B 合法逻辑块和保护组统计。',
         '- CHEMBL45183.json：上述示例的实际 IDs 和内容。',
         '- summary.json / verification.json：各数据划分结果和独立检查。','']
    (folder/'REPORT.zh.md').write_text('\n'.join(md),encoding='utf-8')
    print(json.dumps(dict(audit=audit,all={k:v['mean'] for k,v in summary['all'].items()},
        example_lengths={e['name']:e['length'] for e in examples}),indent=2))


if __name__=='__main__':
    main()
