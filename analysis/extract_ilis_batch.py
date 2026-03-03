#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from collections import OrderedDict
from pathlib import Path

import numpy as np
from scipy.spatial.distance import pdist, squareform


def transform_pae_matrix(pae_matrix: np.ndarray, pae_cutoff: float) -> np.ndarray:
    transformed = np.zeros_like(pae_matrix, dtype=float)
    within = pae_matrix < pae_cutoff
    transformed[within] = 1.0 - (pae_matrix[within] / pae_cutoff)
    return transformed


def calculate_mean_lis(matrix: np.ndarray, subunit_sizes: list[int]) -> np.ndarray:
    cum = np.cumsum(subunit_sizes)
    starts = np.concatenate(([0], cum[:-1]))
    out = np.zeros((len(subunit_sizes), len(subunit_sizes)), dtype=float)

    for i in range(len(subunit_sizes)):
        for j in range(len(subunit_sizes)):
            block = matrix[starts[i] : cum[i], starts[j] : cum[j]]
            vals = block[block > 0]
            out[i, j] = float(vals.mean()) if vals.size else 0.0
    return out


def calculate_contact_map(cif_file: Path, distance_threshold: float = 8.0) -> np.ndarray:
    residue_lines: list[str] = []
    with cif_file.open() as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            if line.startswith("ATOM") and ("CB" in line or ("GLY" in line and "CA" in line)):
                residue_lines.append(line)
            elif line.startswith("ATOM") and "P   " in line:
                residue_lines.append(line)
            elif line.startswith("HETATM"):
                residue_lines.append(line)

    if not residue_lines:
        return np.zeros((0, 0), dtype=int)

    rows: list[list[object]] = []
    for line in residue_lines:
        parts: list[object] = []
        for p in line.split():
            try:
                parts.append(float(p))
            except ValueError:
                parts.append(p)
        rows.append(parts)

    max_len = max(len(r) for r in rows)
    padded = [r + [""] * (max_len - len(r)) for r in rows]
    arr = np.array(padded, dtype=object)

    if arr.shape[1] < 14:
        return np.zeros((arr.shape[0], arr.shape[0]), dtype=int)

    coords = arr[:, 11:14].astype(float)
    distances = squareform(pdist(coords)) if len(coords) > 1 else np.zeros((len(coords), len(coords)))
    atom_col = arr[:, 3].astype(str)
    has_phosphorus = np.array(["P" in a for a in atom_col])
    adjusted = np.where(has_phosphorus[:, None] | has_phosphorus[None, :], distances - 4.0, distances)
    return (adjusted < distance_threshold).astype(int)


def ordered_chain_sizes(token_chain_ids: list[str]) -> list[int]:
    counts: OrderedDict[str, int] = OrderedDict()
    for cid in token_chain_ids:
        counts[cid] = counts.get(cid, 0) + 1
    return list(counts.values())


def pick_file(d: Path, preferred: str, pattern: str) -> Path | None:
    p = d / preferred
    if p.exists():
        return p
    matches = sorted(d.glob(pattern))
    return matches[0] if matches else None


def interchain_scalar(mat: np.ndarray, mode: str = "max") -> float:
    if mat.size == 0 or mat.shape[0] < 2:
        return 0.0
    mask = ~np.eye(mat.shape[0], dtype=bool)
    vals = mat[mask]
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return 0.0
    return float(np.max(vals) if mode == "max" else np.mean(vals))


def compute_ilis_for_dir(run_dir: Path, pae_cutoff: float, distance_cutoff: float) -> float | None:
    conf = pick_file(run_dir, "sp_confidences.json", "*_confidences.json")
    cif = pick_file(run_dir, "sp_model.cif", "*_model.cif")
    if conf is None or cif is None:
        return None

    try:
        with conf.open() as f:
            data = json.load(f)
    except json.JSONDecodeError:
        return None
    except OSError:
        return None

    pae = np.array(data.get("pae", []), dtype=float)
    token_chain_ids = data.get("token_chain_ids", [])
    if pae.size == 0 or not token_chain_ids:
        return None

    subunit_sizes = ordered_chain_sizes(token_chain_ids)
    if sum(subunit_sizes) != pae.shape[0]:
        return None

    transformed = transform_pae_matrix(pae, pae_cutoff)
    mean_lis = calculate_mean_lis(transformed, subunit_sizes)

    contact_map = calculate_contact_map(cif, distance_cutoff)
    if contact_map.shape != transformed.shape:
        n = min(contact_map.shape[0], transformed.shape[0])
        if n == 0:
            return None
        trimmed_t = transformed[:n, :n]
        trimmed_c = contact_map[:n, :n]
        transformed = trimmed_t
        mean_lis = calculate_mean_lis(trimmed_t, subunit_sizes)
        combined = np.where((trimmed_t > 0) & (trimmed_c == 1), trimmed_t, 0.0)
    else:
        combined = np.where((transformed > 0) & (contact_map == 1), transformed, 0.0)

    mean_clis = calculate_mean_lis(combined, subunit_sizes)
    ilis = np.sqrt(np.clip(mean_lis * mean_clis, 0.0, None))
    return interchain_scalar(ilis, mode="max")


def update_score_summary_csv(score_csv: Path, directory: str, ilis_value: float | None) -> None:
    row = {
        "directory": directory,
        "ipTM": "NA",
        "IPSAE": "NA",
        "LIS": "NA",
        "iLIS": "NA" if ilis_value is None else f"{ilis_value:.6f}",
    }

    if score_csv.exists():
        with score_csv.open(newline="") as f:
            reader = csv.DictReader(f)
            existing = next(reader, None)
            if existing:
                for key in ["directory", "ipTM", "IPSAE", "LIS", "iLIS"]:
                    if key in existing and existing[key] != "":
                        row[key] = existing[key]
                if ilis_value is not None:
                    row["iLIS"] = f"{ilis_value:.6f}"

    score_csv.parent.mkdir(parents=True, exist_ok=True)
    with score_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["directory", "ipTM", "IPSAE", "LIS", "iLIS"])
        writer.writeheader()
        writer.writerow(row)


def rebuild_all_scores(final_dirs: Path, output_file: Path) -> None:
    rows: list[dict[str, str]] = []
    for d in sorted(final_dirs.iterdir()):
        if not d.is_dir():
            continue
        f = d / "score_summary.csv"
        if not f.exists():
            continue
        with f.open(newline="") as fh:
            r = csv.DictReader(fh)
            row = next(r, None)
            if not row:
                continue
            rows.append(
                {
                    "directory": row.get("directory", d.name),
                    "ipTM": row.get("ipTM", "NA"),
                    "IPSAE": row.get("IPSAE", "NA"),
                    "LIS": row.get("LIS", "NA"),
                    "iLIS": row.get("iLIS", "NA"),
                }
            )

    with output_file.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["directory", "ipTM", "IPSAE", "LIS", "iLIS"])
        w.writeheader()
        w.writerows(rows)


def run_batch(project_root: Path, pae_cutoff: float, distance_cutoff: float) -> tuple[int, int]:
    out_dirs = project_root / "out_dirs"
    final_dirs = project_root / "final_dirs"
    final_dirs.mkdir(exist_ok=True)

    processed = 0
    missing = 0
    errors: list[str] = []

    for run_dir in sorted([d for d in out_dirs.iterdir() if d.is_dir()]):
        try:
            ilis_value = compute_ilis_for_dir(run_dir, pae_cutoff, distance_cutoff)
        except Exception as e:
            ilis_value = None
            errors.append(f"{run_dir.name}: {type(e).__name__}: {e}")
        if ilis_value is None:
            missing += 1
            # Keep lightweight breadcrumbs for expected data issues.
            conf = pick_file(run_dir, "sp_confidences.json", "*_confidences.json")
            model = pick_file(run_dir, "sp_model.cif", "*_model.cif")
            if conf is None or model is None:
                errors.append(f"{run_dir.name}: missing required files (conf={conf is not None}, model={model is not None})")
        score_csv = final_dirs / run_dir.name / "score_summary.csv"
        update_score_summary_csv(score_csv, run_dir.name, ilis_value)
        processed += 1

    rebuild_all_scores(final_dirs, final_dirs / "all_scores.csv")
    if errors:
        with (final_dirs / "_ilis_errors.log").open("w") as f:
            for e in errors:
                f.write(e + "\n")
    return processed, missing


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Batch compute iLIS for out_dirs and update final_dirs score CSVs.")
    parser.add_argument(
        "--project-root",
        type=Path,
        default=Path("/Users/caysonhamilton/Classes/BIO465/AlphaFold_Stuff"),
        help="Project root containing out_dirs and final_dirs",
    )
    parser.add_argument("--pae-cutoff", type=float, default=12.0)
    parser.add_argument("--distance-cutoff", type=float, default=8.0)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    processed, missing = run_batch(args.project_root, args.pae_cutoff, args.distance_cutoff)
    print(f"Processed directories: {processed}")
    print(f"Directories without computable iLIS: {missing}")
    print(f"Updated: {args.project_root / 'final_dirs' / 'all_scores.csv'}")


if __name__ == "__main__":
    main()
