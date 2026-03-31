#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
from pathlib import Path

import numpy as np
import pandas as pd


METRICS = ["ipTM", "IPSAE", "LIS", "iLIS"]


def normalize_weights(raw_weights: str) -> dict[str, float]:
    parts = [p.strip() for p in raw_weights.split(",")]
    if len(parts) != len(METRICS):
        raise ValueError(
            f"--weights must have {len(METRICS)} values in this order: {', '.join(METRICS)}"
        )
    values = np.array([float(x) for x in parts], dtype=float)
    if np.any(values < 0):
        raise ValueError("Weights must be non-negative.")
    if values.sum() == 0:
        raise ValueError("At least one weight must be > 0.")
    values = values / values.sum()
    return dict(zip(METRICS, values.tolist()))


def rank_scores(df: pd.DataFrame, weights: dict[str, float]) -> pd.DataFrame:
    ranked = df.copy()

    for metric in METRICS:
        ranked[metric] = pd.to_numeric(ranked[metric], errors="coerce")

    ranked = ranked.dropna(subset=["directory"]).copy()

    # Per-metric ranks and percentiles (higher score = better rank/percentile)
    for metric in METRICS:
        ranked[f"{metric}_rank"] = ranked[metric].rank(method="min", ascending=False)
        ranked[f"{metric}_pct"] = ranked[metric].rank(method="average", ascending=True, pct=True)

    rank_cols = [f"{m}_rank" for m in METRICS]
    pct_cols = [f"{m}_pct" for m in METRICS]

    ranked["avg_rank"] = ranked[rank_cols].mean(axis=1, skipna=True)
    ranked["avg_percentile"] = ranked[pct_cols].mean(axis=1, skipna=True)

    # Weighted percentile score (0-1), robust to missing metrics per row
    weighted_num = np.zeros(len(ranked), dtype=float)
    weighted_den = np.zeros(len(ranked), dtype=float)
    for metric in METRICS:
        col = ranked[f"{metric}_pct"].to_numpy(dtype=float)
        w = weights[metric]
        mask = np.isfinite(col)
        weighted_num[mask] += w * col[mask]
        weighted_den[mask] += w
    with np.errstate(invalid="ignore", divide="ignore"):
        ranked["weighted_percentile"] = np.where(weighted_den > 0, weighted_num / weighted_den, np.nan)

    # Weighted z-score composite for additional separation
    weighted_z_num = np.zeros(len(ranked), dtype=float)
    weighted_z_den = np.zeros(len(ranked), dtype=float)
    for metric in METRICS:
        s = ranked[metric]
        mu = s.mean(skipna=True)
        sigma = s.std(skipna=True)
        if pd.isna(sigma) or sigma == 0:
            z = pd.Series(np.nan, index=ranked.index)
        else:
            z = (s - mu) / sigma
        ranked[f"{metric}_z"] = z
        zvals = z.to_numpy(dtype=float)
        w = weights[metric]
        mask = np.isfinite(zvals)
        weighted_z_num[mask] += w * zvals[mask]
        weighted_z_den[mask] += w
    with np.errstate(invalid="ignore", divide="ignore"):
        ranked["weighted_z"] = np.where(weighted_z_den > 0, weighted_z_num / weighted_z_den, np.nan)

    ranked["overall_rank"] = ranked["weighted_percentile"].rank(method="min", ascending=False)

    ranked = ranked.sort_values(
        by=["overall_rank", "weighted_percentile", "avg_percentile", "avg_rank"],
        ascending=[True, False, False, True],
    )

    return ranked


def add_kinase_ids(
    ranked: pd.DataFrame,
    mapping_csv: Path | None,
    uniprot_col: str = "UniprotID",
    primary_label_col: str = "HGNC Name",
    fallback_label_col: str = "xName",
) -> pd.DataFrame:
    out = ranked.copy()
    out["kinase_id"] = np.nan
    if mapping_csv is None or not mapping_csv.exists():
        return out

    m = pd.read_csv(mapping_csv)
    if uniprot_col not in m.columns:
        return out

    if primary_label_col in m.columns:
        primary = m[primary_label_col].astype(str).str.strip()
    else:
        primary = pd.Series([""] * len(m))

    if fallback_label_col in m.columns:
        fallback = m[fallback_label_col].astype(str).str.strip()
    else:
        fallback = pd.Series([""] * len(m))

    label = primary.where(primary != "", fallback)
    label = label.where(label != "", np.nan)

    map_df = pd.DataFrame(
        {
            "directory": m[uniprot_col].astype(str).str.strip(),
            "kinase_id": label,
        }
    ).dropna(subset=["directory"]).drop_duplicates(subset=["directory"], keep="first")

    out = out.merge(map_df, on="directory", how="left", suffixes=("", "_mapped"))
    if "kinase_id_mapped" in out.columns:
        out["kinase_id"] = out["kinase_id_mapped"]
        out = out.drop(columns=["kinase_id_mapped"])
    return out


def write_markdown_report(ranked: pd.DataFrame, out_path: Path, top_n: int, weights: dict[str, float]) -> None:
    lines: list[str] = []
    lines.append("# Score Ranking Report")
    lines.append("")
    lines.append("## Composite Ranking Method")
    lines.append("")
    lines.append("- Composite: `weighted_percentile` (higher is better)")
    lines.append("- Secondary: `avg_percentile`, `avg_rank`")
    lines.append(
        "- Weights (`ipTM, IPSAE, LIS, iLIS`): "
        + ", ".join(f"{m}={weights[m]:.3f}" for m in METRICS)
    )
    lines.append("")

    top = ranked[
        ["overall_rank", "directory", "kinase_id", "weighted_percentile", "avg_percentile", "avg_rank"] + METRICS
    ].head(top_n)
    lines.append(f"## Top {top_n} Overall")
    lines.append("")
    lines.append(top.to_markdown(index=False, floatfmt=".6f"))
    lines.append("")

    for metric in METRICS:
        mtop = (
            ranked[["directory", "kinase_id", metric, f"{metric}_rank", f"{metric}_pct"]]
            .dropna(subset=[metric])
            .sort_values(metric, ascending=False)
            .head(top_n)
        )
        lines.append(f"## Top {top_n} by {metric}")
        lines.append("")
        lines.append(mtop.to_markdown(index=False, floatfmt=".6f"))
        lines.append("")

    out_path.write_text("\n".join(lines))


def _fmt(x: object) -> str:
    if isinstance(x, (int, float, np.floating)):
        if pd.isna(x):
            return ""
        return f"{float(x):.4f}"
    return str(x)


def _df_to_html_table(df: pd.DataFrame) -> str:
    header = "".join(f"<th>{c}</th>" for c in df.columns)
    rows = []
    for _, row in df.iterrows():
        rendered_cells = []
        for col, val in zip(df.columns.tolist(), row.tolist()):
            if col == "PAE":
                rendered_cells.append(f"<td>{val}</td>")
            else:
                rendered_cells.append(f"<td>{_fmt(val)}</td>")
        cells = "".join(rendered_cells)
        rows.append(f"<tr>{cells}</tr>")
    body = "".join(rows)
    return f"<table><thead><tr>{header}</tr></thead><tbody>{body}</tbody></table>"


def _with_pae_thumbnail(df: pd.DataFrame, final_dirs_root: Path) -> pd.DataFrame:
    out = df.copy()
    if "directory" not in out.columns:
        return out

    def _thumb_html(directory: str) -> str:
        img_path = final_dirs_root / str(directory) / "pae.png"
        if not img_path.exists():
            return ""
        src = os.path.relpath(img_path, start=(final_dirs_root / "rankings")).replace("\\", "/")
        return (
            f'<a href="{src}" target="_blank" rel="noopener noreferrer">'
            f'<img src="{src}" alt="pae" style="width:90px;height:auto;border:1px solid #ccc;" />'
            f"</a>"
        )

    if "kinase_id" in out.columns:
        insert_at = out.columns.get_loc("kinase_id") + 1
    else:
        insert_at = out.columns.get_loc("directory") + 1
    out.insert(insert_at, "PAE", out["directory"].map(_thumb_html))
    return out


def write_html_report(
    ranked: pd.DataFrame, out_path: Path, top_n: int, weights: dict[str, float], final_dirs_root: Path
) -> None:
    top = ranked[
        ["overall_rank", "directory", "kinase_id", "weighted_percentile", "avg_percentile", "avg_rank"] + METRICS
    ].head(top_n)
    top = _with_pae_thumbnail(top, final_dirs_root)

    sections = [
        "<h1>Score Ranking Report</h1>",
        "<p><b>Composite:</b> weighted_percentile (higher is better)</p>",
        "<p><b>Weights:</b> "
        + ", ".join(f"{m}={weights[m]:.3f}" for m in METRICS)
        + "</p>",
        f"<h2>Top {top_n} Overall</h2>",
        _df_to_html_table(top),
    ]

    for metric in METRICS:
        mtop = (
            ranked[["directory", "kinase_id", metric, f"{metric}_rank", f"{metric}_pct"]]
            .dropna(subset=[metric])
            .sort_values(metric, ascending=False)
            .head(top_n)
        )
        mtop = _with_pae_thumbnail(mtop, final_dirs_root)
        sections.append(f"<h2>Top {top_n} by {metric}</h2>")
        sections.append(_df_to_html_table(mtop))

    html = f"""<!doctype html>
<html>
<head>
  <meta charset="utf-8" />
  <title>Score Ranking Report</title>
  <style>
    body {{ font-family: Arial, sans-serif; margin: 24px; color: #111; }}
    h1, h2 {{ margin-bottom: 8px; }}
    table {{ border-collapse: collapse; width: 100%; margin-bottom: 22px; table-layout: fixed; }}
    th, td {{ border: 1px solid #ddd; padding: 6px 8px; font-size: 12px; text-align: right; white-space: nowrap; }}
    th:first-child, td:first-child {{ text-align: left; }}
    th:nth-child(2), td:nth-child(2) {{ text-align: left; }}
    thead th {{ background: #f2f2f2; position: sticky; top: 0; }}
    @media print {{
      body {{ margin: 10mm; }}
      th, td {{ font-size: 10px; }}
      h2 {{ page-break-before: always; }}
    }}
  </style>
</head>
<body>
{''.join(sections)}
</body>
</html>
"""
    out_path.write_text(html)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Rank IDs across ipTM/IPSAE/LIS/iLIS and generate easy-to-review outputs."
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("/Users/caysonhamilton/Classes/BIO465/AlphaFold_Stuff/final_dirs/all_scores.csv"),
        help="Input all_scores.csv path",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("/Users/caysonhamilton/Classes/BIO465/AlphaFold_Stuff/final_dirs/rankings"),
        help="Directory to write ranking outputs",
    )
    parser.add_argument("--top-n", type=int, default=30, help="Top N rows to include in summary tables")
    parser.add_argument(
        "--weights",
        type=str,
        default="1,1,1,1",
        help="Comma-separated weights in order: ipTM,IPSAE,LIS,iLIS",
    )
    parser.add_argument(
        "--mapping-csv",
        type=Path,
        default=Path("/Users/caysonhamilton/Classes/BIO465/AlphaFold_Stuff/kinases_notkl.csv"),
        help="CSV mapping UniProt directory IDs to kinase IDs",
    )
    args = parser.parse_args()

    df = pd.read_csv(args.input)
    weights = normalize_weights(args.weights)
    ranked = rank_scores(df, weights)
    ranked = add_kinase_ids(ranked, args.mapping_csv)

    args.output_dir.mkdir(parents=True, exist_ok=True)

    full_cols = (
        ["overall_rank", "directory", "kinase_id", "weighted_percentile", "weighted_z", "avg_percentile", "avg_rank"]
        + METRICS
        + [f"{m}_rank" for m in METRICS]
        + [f"{m}_pct" for m in METRICS]
    )
    ranked[full_cols].to_csv(args.output_dir / "all_scores_ranked.csv", index=False)

    ranked[
        ["overall_rank", "directory", "kinase_id", "weighted_percentile", "avg_percentile", "avg_rank"] + METRICS
    ].head(args.top_n).to_csv(args.output_dir / "top_overall.csv", index=False)

    for metric in METRICS:
        ranked[
            ["directory", "kinase_id", metric, f"{metric}_rank", f"{metric}_pct"]
        ].dropna(subset=[metric]).sort_values(metric, ascending=False).head(args.top_n).to_csv(
            args.output_dir / f"top_{metric}.csv", index=False
        )

    write_markdown_report(ranked, args.output_dir / "ranking_report.md", args.top_n, weights)
    write_html_report(
        ranked,
        args.output_dir / "ranking_report.html",
        args.top_n,
        weights,
        args.output_dir.parent,
    )

    with pd.ExcelWriter(args.output_dir / "ranking_report.xlsx") as xw:
        ranked[full_cols].to_excel(xw, sheet_name="all_scores_ranked", index=False)
        ranked[
            ["overall_rank", "directory", "kinase_id", "weighted_percentile", "avg_percentile", "avg_rank"] + METRICS
        ].head(args.top_n).to_excel(xw, sheet_name="top_overall", index=False)
        for metric in METRICS:
            (
                ranked[["directory", "kinase_id", metric, f"{metric}_rank", f"{metric}_pct"]]
                .dropna(subset=[metric])
                .sort_values(metric, ascending=False)
                .head(args.top_n)
                .to_excel(xw, sheet_name=f"top_{metric}", index=False)
            )

    print(f"Wrote: {args.output_dir / 'all_scores_ranked.csv'}")
    print(f"Wrote: {args.output_dir / 'top_overall.csv'}")
    print(f"Wrote: {args.output_dir / 'ranking_report.md'}")
    print(f"Wrote: {args.output_dir / 'ranking_report.html'}")
    print(f"Wrote: {args.output_dir / 'ranking_report.xlsx'}")


if __name__ == "__main__":
    main()
