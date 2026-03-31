#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from plotnine import aes, geom_jitter, geom_violin, ggplot, labs, theme_bw


def build_plot(input_csv: Path, output_png: Path, width: float, height: float, dpi: int) -> None:
    df = pd.read_csv(input_csv)

    score_columns = ["ipTM", "IPSAE", "LIS", "iLIS"]
    missing = [c for c in score_columns if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns in CSV: {missing}")

    plot_df = df[score_columns].copy()
    for col in score_columns:
        plot_df[col] = pd.to_numeric(plot_df[col], errors="coerce")

    long_df = plot_df.melt(var_name="metric", value_name="score").dropna(subset=["score"])
    if long_df.empty:
        raise ValueError("No numeric values found in ipTM/IPSAE/LIS/iLIS after parsing.")

    p = (
        ggplot(long_df, aes(x="metric", y="score", fill="metric"))
        + geom_violin(trim=False, alpha=0.8)
        + geom_jitter(width=0.12, size=1.2, alpha=0.45, color="black")
        + theme_bw()
        + labs(
            title="Score Distributions from all_scores.csv",
            x="Metric",
            y="Score",
        )
    )

    output_png.parent.mkdir(parents=True, exist_ok=True)
    p.save(filename=str(output_png), width=width, height=height, dpi=dpi, verbose=False)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create violin plots for ipTM, IPSAE, LIS, and iLIS from all_scores.csv using plotnine."
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("/Users/caysonhamilton/Classes/BIO465/AlphaFold_Stuff/final_dirs/all_scores.csv"),
        help="Path to all_scores.csv",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("/Users/caysonhamilton/Classes/BIO465/AlphaFold_Stuff/final_dirs/all_scores_violin.png"),
        help="Output PNG path",
    )
    parser.add_argument("--width", type=float, default=8.0, help="Plot width in inches")
    parser.add_argument("--height", type=float, default=5.0, help="Plot height in inches")
    parser.add_argument("--dpi", type=int, default=300, help="Output DPI")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    build_plot(args.input, args.output, args.width, args.height, args.dpi)
    print(f"Saved violin plot to: {args.output}")


if __name__ == "__main__":
    main()
