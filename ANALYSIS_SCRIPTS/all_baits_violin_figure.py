"""
Usage:
    python kinase_violin_plot.py

Each CSV must have columns: directory, ipTM, IPSAE, LIS, iLIS
One violin per file, showing the distribution of per-kinase composite scores
(row-wise mean of the four score columns, ignoring NAs).
"""

import os
import pandas as pd
import plotly.express as px
import plotly.io as pio
pio.renderers.default = "json"  # prevents auto browser open

# change files and labels based on intended usage
FILES = ["smo_all_scores.csv", "gli1_all_scores.csv", "gli2_all_scores.csv", "sufu_all_scores.csv"]

LABELS = {
    "smo":  "SMO C-Term",
    "gli1": "GLI1 N-Term",
    "gli2": "GLI2",
    "sufu": "SUFU",
}


def load_composite_scores(filepath):
    df = pd.read_csv(filepath)
    score_cols = ["ipTM", "IPSAE", "LIS", "iLIS"]
    missing = [c for c in score_cols if c not in df.columns]
    if missing:
        raise ValueError(f"{filepath}: missing columns {missing}")
    df[score_cols] = df[score_cols].apply(pd.to_numeric, errors="coerce")
    composite = df[score_cols].mean(axis=1, skipna=True).dropna()
    key = os.path.basename(filepath).replace("_all_scores.csv", "")
    label = LABELS.get(key, key)
    return label, composite


def lookup_uniprot_names(uniprot_ids):
    import urllib.request, urllib.parse, ssl
    id_str = " OR ".join(f"accession:{uid}" for uid in uniprot_ids)
    params = urllib.parse.urlencode({
        "query": id_str,
        "fields": "accession,gene_names",
        "format": "tsv",
        "size": len(uniprot_ids),
    })
    url = f"https://rest.uniprot.org/uniprotkb/search?{params}"
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE
    name_map = {}
    try:
        with urllib.request.urlopen(url, timeout=15, context=ctx) as resp:
            lines = resp.read().decode().strip().splitlines()
        for line in lines[1:]:
            parts = line.split("\t")
            if len(parts) >= 2:
                accession = parts[0].strip()
                gene = parts[1].strip().split()[0] if parts[1].strip() else accession
                name_map[accession] = gene
    except Exception as e:
        print(f"  Warning: UniProt lookup failed ({e}), falling back to IDs.")
    for uid in uniprot_ids:
        if uid not in name_map:
            name_map[uid] = uid
    return name_map


def plot_violins(files, output_path="kinase_violin_plot.html"):
    rows = []
    ordered_labels = []
    for f in files:
        label, scores = load_composite_scores(f)
        ordered_labels.append(label)
        print(f"  {label}: {len(scores)} kinases, mean={scores.mean():.3f}, median={scores.median():.3f}")
        for val in scores:
            rows.append({"Protein": label, "Composite Score": val})
    df = pd.DataFrame(rows)
    fig = px.violin(
        df, x="Protein", y="Composite Score", color="Protein",
        box=True, points="all",
        category_orders={"Protein": ordered_labels},
        title="Kinase Composite Score Distributions",
        labels={"Composite Score": "Composite Score (mean of ipTM, IPSAE, LIS, iLIS)"},
        color_discrete_sequence=px.colors.qualitative.T10,
    )
    fig.update_traces(pointpos=0, jitter=0.3, marker=dict(size=4, opacity=0.5))
    fig.update_layout(
        showlegend=False,
        yaxis=dict(range=[-0.05, 1.05], gridcolor="lightgrey"),
        xaxis_title=None,
        plot_bgcolor="white",
        font=dict(size=13),
        title_font=dict(size=16, family="Arial Black"),
    )
    fig.write_html(output_path, auto_open=False)
    print(f"Plot saved to: {output_path}")


def plot_top50_kinases(files, output_path="top50_kinases.html"):
    score_cols = ["ipTM", "IPSAE", "LIS", "iLIS"]
    per_bait = {}
    for f in files:
        df = pd.read_csv(f)
        df[score_cols] = df[score_cols].apply(pd.to_numeric, errors="coerce")
        df["composite"] = df[score_cols].mean(axis=1, skipna=True)
        key = os.path.basename(f).replace("_all_scores.csv", "")
        label = LABELS.get(key, key)
        per_bait[label] = df.set_index("directory")["composite"]
    combined = pd.DataFrame(per_bait)
    combined["total"] = combined.sum(axis=1, skipna=True)
    top50 = combined.nlargest(50, "total").drop(columns="total")
    top50 = top50.loc[top50.sum(axis=1).sort_values(ascending=False).index]
    print("  Looking up gene names from UniProt...")
    top50_reset = top50.reset_index().rename(columns={"directory": "Kinase"})
    name_map = lookup_uniprot_names(top50_reset["Kinase"].tolist())
    top50_reset["Kinase"] = top50_reset["Kinase"].map(name_map)
    melted = top50_reset.melt(id_vars="Kinase", var_name="Bait Protein", value_name="Composite Score")
    ordered_baits = [LABELS[os.path.basename(f).replace("_all_scores.csv", "")] for f in files]
    fig = px.bar(
        melted, x="Composite Score", y="Kinase", color="Bait Protein",
        orientation="h",
        title="Top 50 Kinases by Total Composite Score",
        category_orders={"Bait Protein": ordered_baits, "Kinase": top50_reset["Kinase"].tolist()},
        color_discrete_sequence=px.colors.qualitative.T10,
        labels={"Composite Score": "Composite Score (sum across bait proteins)"},
    )
    fig.update_layout(
        barmode="stack", height=1100, plot_bgcolor="white",
        xaxis=dict(gridcolor="lightgrey"), yaxis_title=None,
        font=dict(size=12), title_font=dict(size=16, family="Arial Black"),
        legend=dict(title="Bait Protein"),
    )
    fig.write_html(output_path, auto_open=False)
    print(f"Top-50 plot saved to: {output_path}")


def build_manuscript_figure(files, output_path="manuscript_figure.html"):
    import json
    import plotly.graph_objects as go
    import numpy as np

    score_cols = ["ipTM", "IPSAE", "LIS", "iLIS"]
    per_bait = {}
    violin_rows = []
    ordered_labels = []

    for f in files:
        df = pd.read_csv(f)
        df[score_cols] = df[score_cols].apply(pd.to_numeric, errors="coerce")
        df["composite"] = df[score_cols].mean(axis=1, skipna=True)
        key = os.path.basename(f).replace("_all_scores.csv", "")
        label = LABELS.get(key, key)
        ordered_labels.append(label)
        per_bait[label] = df.set_index("directory")["composite"]
        for val in df["composite"].dropna():
            violin_rows.append({"Protein": label, "Composite Score": round(val, 4)})

    # Top 10
    combined = pd.DataFrame(per_bait)
    combined["Total"] = combined.sum(axis=1, skipna=True)
    top10 = combined.nlargest(10, "Total").copy()
    top10 = top10.sort_values("Total", ascending=False)
    for col in ordered_labels + ["Total"]:
        top10[col] = top10[col].round(3)
    top10 = top10.reset_index().rename(columns={"directory": "Kinase"})

    print("  Looking up gene names from UniProt...")
    uniprot_ids = top10["Kinase"].tolist()
    name_map = lookup_uniprot_names(uniprot_ids)
    top10["Gene"] = top10["Kinase"].map(name_map)

    top10_ids = top10["Kinase"].tolist()
    top10_genes = top10["Gene"].tolist()
    kinase_colors = px.colors.qualitative.D3[:10]
    id_to_color = {uid: kinase_colors[i] for i, uid in enumerate(top10_ids)}
    id_to_gene = dict(zip(top10_ids, top10_genes))

    traces = []
    for bait in ordered_labels:
        bait_scores = per_bait[bait]
        all_vals = bait_scores.dropna().values
        traces.append(go.Violin(
            x=[bait] * len(all_vals),
            y=all_vals,
            name=bait,
            showlegend=False,
            line_color="#aaaaaa",
            fillcolor="rgba(200,200,200,0.4)",
            points="all",
            pointpos=0,
            jitter=0.3,
            marker=dict(color="#cccccc", size=4, opacity=0.5),
            box_visible=False,
            meanline_visible=False,
        ))
        # Sort top10 ascending by score so highest-scoring dots are added last (drawn on top)
        scored_top10 = []
        for uid in top10_ids:
            if uid in bait_scores.index and not pd.isna(bait_scores[uid]):
                scored_top10.append((bait_scores[uid], uid))
        scored_top10.sort(key=lambda x: x[0])  # lowest first, highest last = on top

        for val, uid in scored_top10:
            gene = id_to_gene[uid]
            color = id_to_color[uid]
            traces.append(go.Scatter(
                x=[bait],
                y=[val],
                mode="markers",
                name=gene,
                legendgroup=gene,
                showlegend=False,
                marker=dict(color=color, size=8, line=dict(color="white", width=0.5)),
                hovertemplate=f"<b>{gene}</b><br>{bait}: %{{y:.3f}}<extra></extra>",
            ))

    vfig = go.Figure(data=traces)
    vfig.update_layout(
        violinmode="overlay",
        xaxis=dict(type="category", categoryorder="array", categoryarray=ordered_labels, tickfont=dict(size=18)),
        yaxis=dict(range=[-0.05, 1.05], gridcolor="lightgrey", title="Composite Score", title_font=dict(size=18)),
        xaxis_title=None,
        plot_bgcolor="white",
        font=dict(size=13),
        showlegend=False,
        margin=dict(b=60, r=20),
    )
    violin_json = json.dumps(vfig.to_dict(), default=str)

    header_cells = "".join(f"<th>{c}</th>" for c in ordered_labels + ["Total"])
    table_rows_html = ""
    for i, (_, row) in enumerate(top10.iterrows()):
        color = kinase_colors[i]
        swatch = f"<span style='display:inline-block;width:10px;height:10px;border-radius:50%;background:{color};margin-right:6px;'></span>"
        cells = "".join(f"<td>{row[c]:.3f}</td>" for c in ordered_labels)
        total_cell = f"<td>{row['Total']:.3f}</td>"
        table_rows_html += f"<tr><td>{swatch}{row['Gene']}</td>{cells}{total_cell}</tr>\n"

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8"/>
<title>Manuscript Figure</title>
<script src="https://cdn.plot.ly/plotly-2.32.0.min.js"></script>
<style>
  body {{ font-family: Arial, sans-serif; font-size: 12px; margin: 20px; background: white; color: black; }}
  .layout {{ display: flex; align-items: flex-start; gap: 32px; }}
  #violin-plot {{ width: 680px; flex-shrink: 0; height: 700px; }}
  table {{ border-collapse: collapse; font-size: 12px; }}
  th, td {{ border: 1px solid #999; padding: 5px 10px; text-align: right; }}
  th {{ background: #f0f0f0; font-weight: bold; }}
  td:first-child, th:first-child {{ text-align: left; }}
</style>
</head>
<body>
<div class="layout">
  <div id="violin-plot"></div>
  <table>
    <thead><tr><th>Kinase</th>{header_cells}</tr></thead>
    <tbody>{table_rows_html}</tbody>
  </table>
</div>
<script>
const violinData = {violin_json};
Plotly.newPlot('violin-plot', violinData.data, violinData.layout, {{responsive: false, displayModeBar: false}});
</script>
</body>
</html>"""

    with open(output_path, "w") as fh:
        fh.write(html)
    print(f"Manuscript figure saved to: {output_path}")


if __name__ == "__main__":
    print("Loading files...")
    plot_violins(FILES)
    print("\nBuilding top-50 kinase chart...")
    plot_top50_kinases(FILES)
    print("\nBuilding manuscript figure...")
    build_manuscript_figure(FILES)