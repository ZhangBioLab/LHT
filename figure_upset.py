#!/usr/bin/env python3
# -*- coding: utf-8 -*-



import argparse
import sys
from pathlib import Path
from typing import Dict, Set, List

import pandas as pd
import matplotlib.pyplot as plt


def read_set(file_path: Path, encoding: str = "utf-8") -> Set[str]:
    elems: Set[str] = set()
    with file_path.open("r", encoding=encoding, errors="ignore") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            elems.add(s)
    return elems


def to_membership_matrix(sets_dict: Dict[str, Set[str]]) -> pd.DataFrame:
    """{set_name: set(elements)} -> DataFrame[bool] (index=元素, columns=集合)."""
    if not sets_dict:
        return pd.DataFrame()
    all_elems: List[str] = sorted(set().union(*sets_dict.values()))
    df = pd.DataFrame(index=all_elems, columns=list(sets_dict.keys()), dtype=bool)
    for sname, sset in sets_dict.items():
        df[sname] = [e in sset for e in all_elems]
    return df


def calc_intersection_counts(df_bool: pd.DataFrame) -> pd.DataFrame:
    
    if df_bool.empty:
        return pd.DataFrame(columns=["combination", "size"])
    combos: List[str] = []
    colnames = list(df_bool.columns)
    for _, row in df_bool.iterrows():
        involved = [colnames[i] for i, v in enumerate(row.values) if v]
        if involved:
            combos.append("&".join(involved))
    s = pd.Series(combos).value_counts()
    out = s.rename_axis("combination").reset_index(name="size")
    return out.sort_values("size", ascending=False, ignore_index=True)


def plot_with_upsetplot(df_bool: pd.DataFrame, out_path: Path, dpi: int):
    
    from upsetplot import UpSet, from_indicators

    
    target_top_down = [
        "Sec",
        "SeU",
        "Se-cofactor",
        "Selenoneine",
        "Ovoselenol",
    ]
    # 自下而上应为小->大
    desired_bottom_up = list(reversed(target_top_down))

    
    set_sizes = df_bool.astype(int).sum(axis=0).to_dict()

    
    size_then_rank = []
    rank_map = {name: i for i, name in enumerate(desired_bottom_up)}
    for name in df_bool.columns:
        size = set_sizes.get(name, 0)
        rank = rank_map.get(name, len(rank_map) + 999)
        size_then_rank.append((size, rank, name))
    ordered_cols = [name for _, _, name in sorted(size_then_rank, key=lambda x: (x[0], x[1]))]

    df_ordered = df_bool[ordered_cols]
    series = from_indicators(ordered_cols, df_ordered)

    plt.figure(figsize=(10, 6))
    u = UpSet(
        series,
        show_counts="%d",
        sort_categories_by=None,   
        sort_by="cardinality"     
    )
    u.plot()

    
    fig = plt.gcf()
    for ax in fig.axes:
        if (ax.get_ylabel() or "").strip().lower() in {"intersection size", "intersections"}:
            ax.set_ylabel("Distribution of organisms", labelpad=20)

    
    base = out_path.with_suffix("")
    plt.savefig(base.with_suffix(".png"), dpi=dpi, bbox_inches="tight")
    plt.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
    plt.savefig(base.with_suffix(".svg"), bbox_inches="tight")
    plt.close()

    print("[INFO] 期望自上到下：", target_top_down)
    print("[INFO] 实际自上到下：", list(reversed(ordered_cols)))
    print(f"[INFO] 已输出: {base.with_suffix('.png')}, {base.with_suffix('.pdf')}, {base.with_suffix('.svg')}")


def main():
    parser = argparse.ArgumentParser(description="Make an UpSet plot from multiple set files.")
    parser.add_argument("--dir", required=True, help="包含集合 txt 文件的目录路径")
    parser.add_argument(
        "--files",
        nargs="+",
        default=[
            "Sec.txt",
            "SeU.txt",
            "Se-cofactor.txt",
            "Selenoneine.txt",
            "Ovoselenol.txt",
        ],
        help="集合文件名列表",
    )
    parser.add_argument("--out", default="upset_plot.png", help="输出文件名前缀（会导出 png/pdf/svg）")
    parser.add_argument("--matrix_csv", default="membership_matrix.csv", help="导出的 0/1 矩阵 CSV 文件名")
    parser.add_argument("--intersections_csv", default="intersection_counts.csv", help="导出的交集计数 CSV 文件名")
    parser.add_argument("--dpi", type=int, default=200, help="PNG 图片 DPI")
    parser.add_argument("--encoding", default="utf-8", help="文本文件编码")
    args = parser.parse_args()

    in_dir = Path(args.dir)
    out_path = Path(args.out)

    
    sets_dict: Dict[str, Set[str]] = {}
    for fname in args.files:
        fpath = in_dir / fname
        if not fpath.exists():
            print(f"[ERROR] 文件不存在: {fpath}", file=sys.stderr)
            sys.exit(1)
        set_name = Path(fname).stem
        elems = read_set(fpath, encoding=args.encoding)
        print(f"[INFO] 读取 {set_name}: {len(elems)} 个元素")
        sets_dict[set_name] = elems

    
    df_bool = to_membership_matrix(sets_dict)
    df_bool.astype(int).to_csv(args.matrix_csv)
    print(f"[INFO] 已导出元素×集合矩阵: {Path(args.matrix_csv).resolve()}")

   
    inter_df = calc_intersection_counts(df_bool)
    inter_df.to_csv(args.intersections_csv, index=False)
    print(f"[INFO] 已导出交集计数: {Path(args.intersections_csv).resolve()}")

    
    try:
        import upsetplot  # noqa: F401
    except ImportError:
        print("[ERROR] 需要安装 'upsetplot'：pip install upsetplot")
        sys.exit(1)

    plot_with_upsetplot(df_bool, out_path, args.dpi)


if __name__ == "__main__":
    main()

