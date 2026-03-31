#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path


def pick_file(run_dir: Path, preferred: str, pattern: str) -> Path | None:
    p = run_dir / preferred
    if p.exists():
        return p
    matches = sorted(run_dir.glob(pattern))
    return matches[0] if matches else None


class ProgressLogger:
    def __init__(self, log_file: Path):
        self.log_file = log_file
        self.log_file.parent.mkdir(parents=True, exist_ok=True)
        with self.log_file.open("w") as f:
            f.write("")

    def _ts(self) -> str:
        return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    def log(self, message: str) -> None:
        line = f"[{self._ts()}] {message}"
        print(line)
        with self.log_file.open("a") as f:
            f.write(line + "\n")


def extract_iptm(summary_file: Path) -> str:
    try:
        with summary_file.open() as f:
            data = json.load(f)
    except Exception:
        return "NA"

    if isinstance(data.get("iptm"), (int, float)):
        return str(data["iptm"])

    mat = data.get("chain_pair_iptm")
    if isinstance(mat, list):
        vals = []
        for i, row in enumerate(mat):
            if not isinstance(row, list):
                continue
            for j, v in enumerate(row):
                if i != j and isinstance(v, (int, float)):
                    vals.append(float(v))
        if vals:
            return str(sum(vals) / len(vals))
    return "NA"


def run_ipsae(run_dir: Path, ipsae_script: Path, conf_file: Path, model_file: Path, pae_cutoff: int, dist_cutoff: int) -> tuple[str, str, str | None]:
    try:
        subprocess.run(
            [sys.executable, str(ipsae_script), conf_file.name, model_file.name, str(pae_cutoff), str(dist_cutoff)],
            cwd=run_dir,
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
    except subprocess.CalledProcessError:
        return "NA", "NA", "ipsae.py failed"

    txt_file = run_dir / f"{model_file.stem}_{pae_cutoff:02d}_{dist_cutoff:02d}.txt"
    if not txt_file.exists():
        return "NA", "NA", f"{txt_file.name} missing"

    rows: list[tuple[str, float, float]] = []
    with txt_file.open() as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("Chn1"):
                continue
            parts = s.split()
            if len(parts) < 13:
                continue
            score_type = parts[4]
            if score_type not in ("max", "asym"):
                continue
            try:
                ipsae = float(parts[5])
                lis = float(parts[12])
            except ValueError:
                continue
            rows.append((score_type, ipsae, lis))

    if not rows:
        return "NA", "NA", "ipsae output parse failed"

    max_rows = [r for r in rows if r[0] == "max"]
    chosen = max(max_rows, key=lambda x: x[1]) if max_rows else max(rows, key=lambda x: x[1])
    return f"{chosen[1]}", f"{chosen[2]}", None


def write_score_summary(score_file: Path, row: dict[str, str]) -> None:
    score_file.parent.mkdir(parents=True, exist_ok=True)
    with score_file.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["directory", "ipTM", "IPSAE", "LIS", "iLIS"])
        w.writeheader()
        w.writerow(row)


def combine_all_scores(final_dirs: Path) -> None:
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
    out = final_dirs / "all_scores.csv"
    with out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["directory", "ipTM", "IPSAE", "LIS", "iLIS"])
        w.writeheader()
        w.writerows(rows)


def process_out_dirs(project_root: Path, pae_cutoff: int, dist_cutoff: int, logger: ProgressLogger) -> Path:
    out_dirs = project_root / "out_dirs"
    final_dirs = project_root / "final_dirs"
    ipsae_script = project_root / "IPSAE" / "ipsae.py"
    final_dirs.mkdir(exist_ok=True)

    run_dirs = sorted([d for d in out_dirs.iterdir() if d.is_dir()])
    total = len(run_dirs)
    errors: list[str] = []
    logger.log(f"Stage process_out_dirs: starting {total} directories")

    for idx, run_dir in enumerate(run_dirs, start=1):
        pct = (idx / total) * 100 if total else 100.0
        logger.log(f"[{idx}/{total} | {pct:5.1f}%] Processing {run_dir.name}")
        target = final_dirs / run_dir.name
        target.mkdir(parents=True, exist_ok=True)

        conf = pick_file(run_dir, "sp_confidences.json", "*_confidences.json")
        summ = pick_file(run_dir, "sp_summary_confidences.json", "*_summary_confidences.json")
        model = pick_file(run_dir, "sp_model.cif", "*_model.cif")

        iptm, ipsae, lis = "NA", "NA", "NA"
        ilis = "NA"

        if model is not None:
            shutil.copy2(model, target / "model.cif")
        else:
            logger.log(f"{run_dir.name}: model file missing")

        # Assume pae.png is already present in out_dirs/<ID>/ (generated remotely).
        pae_src = run_dir / "pae.png"
        if pae_src.exists():
            shutil.copy2(pae_src, target / "pae.png")
        else:
            errors.append(f"{run_dir.name}: pae.png missing")
            logger.log(f"{run_dir.name}: pae.png missing")

        if conf is not None and model is not None:
            ipsae, lis, err2 = run_ipsae(run_dir, ipsae_script, conf, model, pae_cutoff, dist_cutoff)
            if err2:
                errors.append(f"{run_dir.name}: {err2}")
                logger.log(f"{run_dir.name}: {err2}")
        else:
            errors.append(f"{run_dir.name}: missing files conf/model")
            logger.log(f"{run_dir.name}: missing files conf/model")

        if summ is not None:
            iptm = extract_iptm(summ)
        else:
            errors.append(f"{run_dir.name}: summary_confidences missing")
            logger.log(f"{run_dir.name}: summary_confidences missing")

        write_score_summary(
            target / "score_summary.csv",
            {"directory": run_dir.name, "ipTM": iptm, "IPSAE": ipsae, "LIS": lis, "iLIS": ilis},
        )

    combine_all_scores(final_dirs)
    logger.log("Stage process_out_dirs: rebuilt final_dirs/all_scores.csv")

    log_file = final_dirs / "_pipeline_process_log.txt"
    with log_file.open("w") as f:
        for e in errors:
            f.write(e + "\n")
    logger.log(f"Stage process_out_dirs: complete with {len(errors)} warnings/errors")
    return log_file


def run_cmd(cmd: list[str], cwd: Path, logger: ProgressLogger, label: str) -> None:
    logger.log(f"Stage {label}: start -> {' '.join(cmd)}")
    subprocess.run(cmd, cwd=cwd, check=True)
    logger.log(f"Stage {label}: complete")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="End-to-end scoring pipeline for AlphaFold_Stuff-style directory structure."
    )
    parser.add_argument(
        "--project-root",
        type=Path,
        default=Path("/Users/caysonhamilton/Classes/BIO465/AlphaFold_Stuff"),
        help="Root containing out_dirs/, final_dirs/, IPSAE/, analysis/",
    )
    parser.add_argument("--pae-cutoff", type=int, default=10, help="PAE cutoff for ipsae.py")
    parser.add_argument("--dist-cutoff", type=int, default=15, help="Distance cutoff for ipsae.py")
    parser.add_argument("--skip-process", action="store_true", help="Skip out_dirs -> final_dirs base processing")
    parser.add_argument("--skip-ilis", action="store_true", help="Skip iLIS extraction/update")
    parser.add_argument("--skip-plots", action="store_true", help="Skip plot generation")
    parser.add_argument("--skip-rankings", action="store_true", help="Skip ranking/report generation")
    parser.add_argument(
        "--log-file",
        type=Path,
        default=None,
        help="Optional log file path (default: final_dirs/logs/pipeline_<timestamp>.log)",
    )
    args = parser.parse_args()

    analysis_dir = args.project_root / "analysis"
    final_dirs = args.project_root / "final_dirs"
    logs_dir = final_dirs / "logs"
    default_log = logs_dir / f"pipeline_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
    logger = ProgressLogger(args.log_file if args.log_file else default_log)
    logger.log("Pipeline start")
    logger.log(f"Project root: {args.project_root}")

    if not args.skip_process:
        log_file = process_out_dirs(args.project_root, args.pae_cutoff, args.dist_cutoff, logger)
        logger.log(f"Processed out_dirs. Error log: {log_file}")

    if not args.skip_ilis:
        run_cmd(
            [sys.executable, str(analysis_dir / "extract_ilis_batch.py"), "--project-root", str(args.project_root)],
            args.project_root,
            logger,
            "iLIS extraction",
        )
        logger.log("Updated iLIS in score_summary.csv files and rebuilt all_scores.csv")

    if not args.skip_plots:
        run_cmd(
            [
                sys.executable,
                str(analysis_dir / "plot_all_scores_violin.py"),
                "--input",
                str(final_dirs / "all_scores.csv"),
                "--output",
                str(final_dirs / "all_scores_violin.png"),
            ],
            args.project_root,
            logger,
            "plot violin",
        )
        run_cmd(
            [
                sys.executable,
                str(analysis_dir / "plot_all_scores_interactive.py"),
                "--input",
                str(final_dirs / "all_scores.csv"),
                "--output",
                str(final_dirs / "all_scores_interactive.html"),
            ],
            args.project_root,
            logger,
            "plot interactive",
        )
        logger.log(f"Plots updated in {final_dirs}")

    if not args.skip_rankings:
        run_cmd(
            [
                sys.executable,
                str(analysis_dir / "rank_score_ids.py"),
                "--input",
                str(final_dirs / "all_scores.csv"),
                "--output-dir",
                str(final_dirs / "rankings"),
                "--mapping-csv",
                str(args.project_root / "kinases_notkl.csv"),
            ],
            args.project_root,
            logger,
            "ranking",
        )
        logger.log(f"Rankings updated in {final_dirs / 'rankings'}")

    logger.log("Pipeline complete.")
    logger.log(f"Main progress log: {logger.log_file}")


if __name__ == "__main__":
    main()
