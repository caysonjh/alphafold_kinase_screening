#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import plotly.express as px
from plotly.subplots import make_subplots


def load_scores(input_csv: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    df = pd.read_csv(input_csv)
    required = ["directory", "ipTM", "IPSAE", "LIS", "iLIS"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    wide = df[["directory", "ipTM", "IPSAE", "LIS", "iLIS"]].copy()
    for col in ["ipTM", "IPSAE", "LIS", "iLIS"]:
        wide[col] = pd.to_numeric(wide[col], errors="coerce")
    wide = wide.dropna(subset=["ipTM", "IPSAE", "LIS", "iLIS"]).reset_index(drop=True)

    long = wide.melt(id_vars="directory", var_name="metric", value_name="score")
    long = long.merge(wide, on="directory", how="left")
    return wide, long


def build_violin_with_points(long_df: pd.DataFrame):
    hover_df = long_df.copy()
    fig = px.violin(
        hover_df,
        x="metric",
        y="score",
        color="metric",
        points="all",
        hover_data={
            "directory": True,
            "metric": True,
            "score": ":.4f",
            "ipTM": ":.4f",
            "IPSAE": ":.4f",
            "LIS": ":.4f",
            "iLIS": ":.4f",
        },
    )
    fig.update_traces(jitter=0.14, marker=dict(size=5, opacity=0.65, line=dict(width=0)))
    fig.update_layout(
        title="Violin Plot with Individual Points (hover for directory)",
        xaxis_title="Metric",
        yaxis_title="Score",
        template="plotly_white",
    )
    return fig


def build_scatter_matrix(wide_df: pd.DataFrame):
    fig = px.scatter_matrix(
        wide_df,
        dimensions=["ipTM", "IPSAE", "LIS", "iLIS"],
        hover_name="directory",
        title="Scatter Matrix (hover to track the same directory across metrics)",
    )
    fig.update_traces(diagonal_visible=False, marker=dict(size=6, opacity=0.55))
    fig.update_layout(template="plotly_white")
    return fig


def write_combined_html(violin_fig, matrix_fig, output_html: Path) -> None:
    output_html.parent.mkdir(parents=True, exist_ok=True)

    # Use a simple HTML wrapper so both figures are in one file.
    violin_html = violin_fig.to_html(full_html=False, include_plotlyjs="cdn")
    matrix_html = matrix_fig.to_html(full_html=False, include_plotlyjs=False)

    html = f"""<!doctype html>
<html lang=\"en\">
<head>
  <meta charset=\"utf-8\" />
  <meta name=\"viewport\" content=\"width=device-width, initial-scale=1\" />
  <title>Interactive all_scores Plots</title>
  <style>
    body {{ font-family: Arial, sans-serif; margin: 24px; background: #f7f7f7; }}
    h1 {{ margin-bottom: 8px; }}
    p {{ margin-top: 0; color: #444; }}
    .panel {{ background: #fff; border: 1px solid #ddd; border-radius: 10px; padding: 12px; margin-bottom: 18px; }}
  </style>
</head>
<body>
  <h1>Interactive Score Visualizations</h1>
  <p>Hover over points to see directory IDs.</p>
  <div class=\"panel\">{violin_html}</div>
  <div class=\"panel\">{matrix_html}</div>
</body>
</html>
"""

    output_html.write_text(html)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create interactive violin + scatter-matrix plots from final_dirs/all_scores.csv"
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
        default=Path("/Users/caysonhamilton/Classes/BIO465/AlphaFold_Stuff/final_dirs/all_scores_interactive.html"),
        help="Output HTML path",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    wide, long = load_scores(args.input)
    violin_fig = build_violin_with_points(long)
    matrix_fig = build_scatter_matrix(wide)
    write_combined_html(violin_fig, matrix_fig, args.output)
    print(f"Saved interactive HTML to: {args.output}")


if __name__ == "__main__":
    main()
