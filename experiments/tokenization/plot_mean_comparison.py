"""Plot saved mean lengths; DemoDiff counts nodes with connectivity kept separately."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties


def draw(data, output, chinese_font=None):
    chinese = chinese_font is not None
    font = FontProperties(fname=str(chinese_font)) if chinese else FontProperties(family="DejaVu Sans")
    labels = [r["label_zh" if chinese else "label"] for r in data["methods"]]
    values = [r["mean"] for r in data["methods"]]
    fig, ax = plt.subplots(figsize=(11, 6.7 + max(0, len(values) - 7) * .5))
    fig.subplots_adjust(left=.37, right=.96, top=.84, bottom=.22)
    bars = ax.barh(range(len(values)), values, height=.57, color=[r["color"] for r in data["methods"]])
    ax.invert_yaxis()
    ax.set_yticks(range(len(labels)), labels, fontproperties=font, fontsize=11)
    ax.set_xlim(0, 69)
    ax.set_xlabel("平均 token 长度" if chinese else "Mean token length", fontproperties=font, fontsize=11, labelpad=10)
    ax.xaxis.grid(True, color="#E7ECEA")
    ax.set_axisbelow(True)
    for spine in ["top", "right", "left"]:
        ax.spines[spine].set_visible(False)
    ax.spines["bottom"].set_color("#CBD4D0")
    ax.tick_params(axis="y", length=0, pad=10)
    for bar, value in zip(bars, values):
        ax.text(value + .8, bar.get_y() + bar.get_height()/2, f"{value:.3f}", va="center", fontsize=11, fontweight="bold", color="#20343C")
    title = "平均 token 长度对比" if chinese else "Mean token lengths"
    subtitle = "同一批 44,329 条边界一致记录" if chinese else "Same 44,329 boundary-preserving records"
    fig.text(.04, .945, title, fontproperties=font, fontsize=20, color="#20343C")
    fig.text(.04, .895, subtitle, fontproperties=font, fontsize=11, color="#52665E")
    note = ("* DemoDiff 仅计节点 token；连接信息由额外边表及连接位点信息保存。\n"
            "其余方法计入文本连接符号。此图不比较总存储量或计算成本。" if chinese else
            "* DemoDiff counts node tokens only; connectivity is stored in extra edge/attachment tables.\n"
            "Other methods count connectivity syntax in the text. This is not a total-storage or compute comparison.")
    fig.text(.04, .055, note, fontproperties=font, fontsize=9, linespacing=1.7, color="#52665E")
    fig.savefig(output, dpi=220, metadata={"Creator": "Molecular tokenizer comparison"})
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--chinese-font", type=Path)
    parser.add_argument("--structural", type=Path, help="Current A/B export; replaces older scheme rows")
    parser.add_argument("--save-data", type=Path, help="Save the exact data plotted")
    args = parser.parse_args()
    data = json.loads(args.data.read_text(encoding="utf-8"))
    if args.structural:
        current = json.loads(args.structural.read_text(encoding="utf-8"))
        data['methods'] = [r for r in data['methods'] if r['key'] in ('atomwise', 'demodiff', 'brics_safe', 'whole_16000', 'whole_32000')]
        for budget in ('16000', '32000'):
            for mode, label, zh, color in (
                ('A', 'Complete syntax blocks', '完整语法块', '#24765B'),
                ('B', 'Protected structures', '受保护子结构', '#68699C'),
            ):
                row = current['budgets'][budget][mode]['all']
                data['methods'].append(dict(key=f'{mode}_{budget}', label=f'{mode} · {label} · {int(budget)//1000}k',
                    label_zh=f'{mode} · {zh} · {int(budget)//1000}k', color=color,
                    records=row['n'], total=row['total'], mean=row['mean'],
                    connectivity='Included in the sequence, together with fallback block markers'))
        data['sources']['structural'] = 'structural-comparison.json'
        data['scope'] = 'Same 44,329 boundary-preserving records; includes TRAIN, VALIDATION and TEST. Current A/B and reference tokenizers.'
    order = ['atomwise', 'demodiff', 'brics_safe', 'whole_16000', 'A_16000', 'B_16000', 'whole_32000', 'A_32000', 'B_32000']
    data['methods'].sort(key=lambda r: order.index(r['key']))
    if args.save_data:
        args.save_data.write_text(json.dumps(data, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    assert len({r['key'] for r in data['methods']}) == len(data['methods'])
    for row in data["methods"]:
        assert row["records"] == data["records"] == 44329
        assert abs(row["mean"] - row["total"] / row["records"]) < 1e-10
    draw(data, args.output, args.chinese_font)


if __name__ == "__main__":
    main()
