import argparse
import hashlib
import json
import sys
from pathlib import Path


def load_tokenizer(project_root, tokenizer_path):
    """加载本地项目的 SAFETokenizer，并统一长度统计配置。"""
    project_root = Path(project_root).resolve()
    tokenizer_path = Path(tokenizer_path).resolve()

    if tokenizer_path.is_dir():
        tokenizer_path = tokenizer_path / "tokenizer.json"

    if not (project_root / "safe" / "tokenizer.py").is_file():
        raise FileNotFoundError(
            f"找不到项目中的 safe/tokenizer.py：{project_root}"
        )
    if not tokenizer_path.is_file():
        raise FileNotFoundError(f"找不到 tokenizer：{tokenizer_path}")

    # 项目根目录用于导入 safe 包。
    # safe 目录兼容当前源码中部分模块使用的直接导入方式。
    sys.path.insert(0, str(project_root / "safe"))
    sys.path.insert(0, str(project_root))

    from safe.tokenizer import SAFETokenizer
    from tokenizers.models import BPE

    # 项目自己的加载方法会恢复自定义 SAFE pre-tokenizer。
    tokenizer = SAFETokenizer.load(str(tokenizer_path))

    tokenizer.tokenizer.no_padding()
    tokenizer.tokenizer.no_truncation()

    # 避免 BPE dropout 使相同文本的统计结果出现随机变化。
    if isinstance(tokenizer.tokenizer.model, BPE):
        tokenizer.tokenizer.model.dropout = None

    return tokenizer, tokenizer_path


def analyze_safe_file(file_path, tokenizer, add_special_tokens=False):
    """流式统计：一行是一条 SAFE 编码，不保存整个语料或长度列表。"""
    file_path = Path(file_path).resolve()
    if not file_path.is_file():
        raise FileNotFoundError(f"找不到 SAFE 文件：{file_path}")

    count = 0
    blank_lines = 0
    total_tokens = 0
    min_tokens = None
    max_tokens = None
    total_unk_tokens = 0
    sequences_with_unk = 0

    unk_id = tokenizer.unk_token_id

    # utf-8-sig 同时兼容普通 UTF-8 和带 BOM 的文件。
    with file_path.open("r", encoding="utf-8-sig") as fin:
        for line_number, line in enumerate(fin, start=1):
            # 只去掉文件行结束符，不改变 SAFE 文本的其他内容。
            safe_text = line.rstrip("\r\n")

            if not safe_text.strip():
                blank_lines += 1
                continue

            try:
                token_ids = tokenizer.encode(
                    safe_text,
                    ids_only=True,
                    add_special_tokens=add_special_tokens,
                )
            except Exception as exc:
                # 不悄悄丢弃失败行，以免统计样本发生变化。
                raise RuntimeError(
                    f"{file_path} 第 {line_number} 行 tokenization 失败"
                ) from exc

            length = len(token_ids)
            if length == 0:
                raise ValueError(
                    f"{file_path} 第 {line_number} 行非空，但 token 数为 0"
                )

            count += 1
            total_tokens += length
            min_tokens = (
                length if min_tokens is None else min(min_tokens, length)
            )
            max_tokens = (
                length if max_tokens is None else max(max_tokens, length)
            )

            if unk_id is not None:
                unk_count = token_ids.count(unk_id)
                total_unk_tokens += unk_count
                sequences_with_unk += int(unk_count > 0)

    if count == 0:
        raise ValueError(f"文件没有可统计的非空 SAFE 编码：{file_path}")

    return {
        "file": str(file_path),
        "num_sequences": count,
        "blank_lines_skipped": blank_lines,
        "total_tokens": total_tokens,
        "mean_tokens": total_tokens / count,
        "min_tokens": min_tokens,
        "max_tokens": max_tokens,
        "unk_tokens": total_unk_tokens,
        "sequences_with_unk": sequences_with_unk,
    }


def main():
    parser = argparse.ArgumentParser(
        description="比较 GraphBPE 与 BRICS 的 SAFE 平均 token 长度"
    )
    parser.add_argument(
        "--project-root",
        type=Path,
        default=Path(r"E:\GraphBpe+SAFE\SAFE-main"),
        help="包含 safe 文件夹的项目根目录",
    )
    parser.add_argument(
        "--tokenizer",
        type=Path,
        required=True,
        help="已训练的 tokenizer JSON 文件，或包含 tokenizer.json 的目录",
    )
    parser.add_argument(
        "--graphbpe", type=Path, required=True,
        help="GraphBPE 生成的 safe.txt",
    )
    parser.add_argument(
        "--brics", type=Path, required=True,
        help="BRICS 生成的 safe.txt",
    )
    parser.add_argument(
        "--add-special-tokens",
        action="store_true",
        help="计入 tokenizer 配置的特殊 token；默认不计入",
    )
    args = parser.parse_args()

    tokenizer, tokenizer_path = load_tokenizer(
        args.project_root, args.tokenizer
    )

    graphbpe = analyze_safe_file(
        args.graphbpe, tokenizer, args.add_special_tokens
    )
    brics = analyze_safe_file(
        args.brics, tokenizer, args.add_special_tokens
    )

    graphbpe_mean = graphbpe["mean_tokens"]
    brics_mean = brics["mean_tokens"]

    report = {
        "tokenizer": str(tokenizer_path),
        "tokenizer_sha256": hashlib.sha256(
            tokenizer_path.read_bytes()
        ).hexdigest(),
        "add_special_tokens": args.add_special_tokens,
        "padding": False,
        "truncation": False,
        "graphbpe": graphbpe,
        "brics": brics,
        "comparison": {
            # 负数表示 GraphBPE 平均序列更短。
            "graphbpe_minus_brics": graphbpe_mean - brics_mean,

            # 正数表示 GraphBPE 相对 BRICS 的平均长度缩短比例。
            "graphbpe_reduction_percent": (
                (brics_mean - graphbpe_mean) / brics_mean * 100
            ),

            # 数量相同也不能证明两份文件包含相同的分子。
            "same_sequence_count": (
                graphbpe["num_sequences"] == brics["num_sequences"]
            ),
            "same_molecule_set_verified": False,
        },
    }

    # 自动创建输出目录
    output_dir = args.project_root.resolve() / "safe" / "data" / "comparison"
    output_dir.mkdir(parents=True, exist_ok=True)

    # 将统计结果保存为 JSON
    output_file = output_dir / "results.json"
    result_json = json.dumps(report, ensure_ascii=False, indent=2)
    output_file.write_text(result_json + "\n", encoding="utf-8")

    # 同时在终端显示结果和保存位置
    print(result_json)
    print(f"\n结果已保存到：{output_file}")

if __name__ == "__main__":
    main()