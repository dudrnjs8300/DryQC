#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import csv
import gzip
import argparse
import subprocess
import tempfile
import hashlib
from collections import defaultdict

ENGINE_REVISION = "FASTQ_V4_ILLUMINA_ALL_TARGETS"

def die(msg):
    print(msg, file=sys.stderr)
    sys.exit(1)

def open_maybe_gz(path, mode="rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)

def fastq_to_fasta(in_path, out_handle, prefix=None):
    n = 0
    with open_maybe_gz(in_path, "rt") as f:
        while True:
            h = f.readline()
            if not h:
                break
            s = f.readline(); p = f.readline(); q = f.readline()
            if not q:
                break
            h = h.strip(); s = s.strip()
            if not h.startswith("@"):
                continue
            rid = h[1:].split()[0]
            if prefix:
                rid = f"{prefix}:{rid}"
            out_handle.write(f">{rid}\n{s}\n")
            n += 1
    return n

def file_sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def ensure_diamond_db(panel_faa):
    dmnd = panel_faa + ".dmnd"
    source_hash_file = dmnd + ".source.sha256"
    panel_hash = file_sha256(panel_faa)
    recorded_hash = ""
    if os.path.exists(source_hash_file):
        with open(source_hash_file, "r", encoding="utf-8") as f:
            recorded_hash = f.read().strip().split()[0] if f.readable() else ""
    if os.path.exists(dmnd) and recorded_hash == panel_hash:
        return dmnd
    tmp_prefix = dmnd + f".tmp.{os.getpid()}"
    tmp_dmnd = tmp_prefix + ".dmnd"
    try:
        subprocess.run(["diamond", "makedb", "--in", panel_faa, "-d", tmp_prefix], check=True)
        os.replace(tmp_dmnd, dmnd)
        tmp_hash = source_hash_file + f".tmp.{os.getpid()}"
        with open(tmp_hash, "w", encoding="utf-8") as f:
            f.write(panel_hash + "  " + os.path.basename(panel_faa) + "\n")
        os.replace(tmp_hash, source_hash_file)
    finally:
        if os.path.exists(tmp_dmnd):
            os.unlink(tmp_dmnd)
    return dmnd

def parse_marker_id(target_id):
    return str(target_id).split("|", 1)[0]

def load_marker_ids_from_panel(panel_faa):
    marker_ids, seen = [], set()
    with open(panel_faa, "r") as f:
        for line in f:
            if line.startswith(">"):
                marker_id = parse_marker_id(line[1:].strip())
                if marker_id and marker_id not in seen:
                    seen.add(marker_id)
                    marker_ids.append(marker_id)
    if not marker_ids:
        die(f"[ERROR] no marker IDs parsed from panel: {panel_faa}")
    return marker_ids

def platform_defaults(platform):
    if platform == "illumina":
        return {"min_pident": 0.80, "min_aln_aa": 25, "min_bitscore": 50.0, "min_reads_per_marker": 5, "low_data_reads": 100000, "max_target_seqs": 0, "threshold_status": "AUDITED_V1_COMPATIBILITY"}
    if platform == "ont":
        return {"min_pident": 0.55, "min_aln_aa": 20, "min_bitscore": 40.0, "min_reads_per_marker": 2, "low_data_reads": 20000, "max_target_seqs": 0, "threshold_status": "AUDITED_V1_COMPATIBILITY"}
    if platform in ("pacbio_hifi", "pacbio_clr"):
        return {"min_pident": 0.55, "min_aln_aa": 20, "min_bitscore": 40.0, "min_reads_per_marker": 2, "low_data_reads": 20000, "max_target_seqs": 0, "threshold_status": "PROVISIONAL_PROFILE_SPECIFIC_CALIBRATION_REQUIRED"}
    if platform == "pacbio":
        return {"min_pident": 0.55, "min_aln_aa": 20, "min_bitscore": 40.0, "min_reads_per_marker": 2, "low_data_reads": 20000, "max_target_seqs": 0, "threshold_status": "LEGACY_GENERIC_PACBIO_DEPRECATED"}
    die(f"[ERROR] unknown platform: {platform}")

def decide_marker_qc_status(detected, required):
    """Marker-recovery classification only."""
    return "PASS" if detected >= required else "FAIL"

def decide_data_sufficiency_status(total_reads, low_data_reads):
    """
    Preserve the historical read-count cutoff as an audit/review flag.
    It does not override marker-recovery QC.
    """
    if total_reads < low_data_reads:
        return "REVIEW_LOW_READ_COUNT_LEGACY"
    return "NOT_FLAGGED_BY_LEGACY_READ_COUNT"

def decide_legacy_qc_status(detected, required, total_reads, low_data_reads):
    """Exact pre-v3 behavior retained for reproducibility."""
    if total_reads < low_data_reads:
        return "LOW_DATA"
    return "PASS" if detected >= required else "FAIL"

def main():
    ap = argparse.ArgumentParser(description="Assembly-free FASTQ -> marker QC via DIAMOND blastx")
    ap.add_argument("--r1", required=True)
    ap.add_argument("--r2", default=None)
    ap.add_argument("--panel", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--platform", choices=["illumina", "ont", "pacbio_hifi", "pacbio_clr", "pacbio"], default="illumina")
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--min-pident", type=float, default=None)
    ap.add_argument("--min-aln-aa", type=int, default=None)
    ap.add_argument("--min-bitscore", type=float, default=None)
    ap.add_argument("--min-reads-per-marker", type=int, default=None)
    ap.add_argument("--required-genes", type=int, default=24)
    ap.add_argument("--max-target-seqs", type=int, default=None,
                    help="DIAMOND target cap; profile default is 0 (all targets) for Illumina and long reads")
    ap.add_argument("--evalue", default="1e-5")
    ap.add_argument("--low-data-reads", type=int, default=None)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    defaults = platform_defaults(args.platform)
    min_pident = args.min_pident if args.min_pident is not None else defaults["min_pident"]
    min_aln_aa = args.min_aln_aa if args.min_aln_aa is not None else defaults["min_aln_aa"]
    min_bitscore = args.min_bitscore if args.min_bitscore is not None else defaults["min_bitscore"]
    min_reads_per_marker = args.min_reads_per_marker if args.min_reads_per_marker is not None else defaults["min_reads_per_marker"]
    low_data_reads = args.low_data_reads if args.low_data_reads is not None else defaults["low_data_reads"]
    max_target_seqs = args.max_target_seqs if args.max_target_seqs is not None else defaults["max_target_seqs"]
    if max_target_seqs < 0:
        die("[ERROR] --max-target-seqs must be >= 0")
    long_read_profile = args.platform in ("ont", "pacbio_hifi", "pacbio_clr", "pacbio")
    if args.platform == "pacbio":
        print("[WARN] generic --platform pacbio is deprecated; specify pacbio_hifi or pacbio_clr", file=sys.stderr)
    if long_read_profile and max_target_seqs > 0:
        print("[WARN] a finite target cap can censor distinct marker families on long reads", file=sys.stderr)

    if not os.path.exists(args.panel): die(f"[ERROR] panel not found: {args.panel}")
    if not os.path.exists(args.r1): die(f"[ERROR] R1 not found: {args.r1}")
    if args.r2 and not os.path.exists(args.r2): die(f"[ERROR] R2 not found: {args.r2}")

    marker_ids = load_marker_ids_from_panel(args.panel)
    ensure_diamond_db(args.panel)

    with tempfile.TemporaryDirectory() as tmpdir:
        fasta_reads = os.path.join(tmpdir, f"{args.sample}.reads.fa")
        total_reads = 0
        with open(fasta_reads, "w") as outfa:
            total_reads += fastq_to_fasta(args.r1, outfa, prefix="R1")
            if args.r2:
                total_reads += fastq_to_fasta(args.r2, outfa, prefix="R2")

        m8 = os.path.join(args.outdir, f"{args.sample}.ribo27.reads.m8")
        cmd = [
            "diamond", "blastx", "--query", fasta_reads, "--db", args.panel + ".dmnd", "--out", m8,
            "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "qstart", "qend", "sstart", "send", "bitscore", "evalue",
            "--max-target-seqs", str(max_target_seqs), "--evalue", str(args.evalue),
            "--threads", str(args.threads), "--more-sensitive"
        ]
        subprocess.run(cmd, check=True)

    detail_csv = os.path.join(args.outdir, f"{args.sample}.ribo27.fastq.detail.csv")
    summary_csv = os.path.join(args.outdir, f"{args.sample}.ribo27.fastq.summary.csv")

    marker_to_reads = defaultdict(set)
    marker_best = {}
    query_report_counts = defaultdict(int)
    reported_hsps = 0

    with open(m8) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 10:
                continue
            qseqid, sseqid, pident, alnlen, qstart, qend, sstart, send, bitscore, evalue = parts
            reported_hsps += 1
            query_report_counts[qseqid] += 1
            pident = float(pident) / 100.0
            alnlen = int(alnlen)
            bitscore = float(bitscore)

            if pident < min_pident or alnlen < min_aln_aa or bitscore < min_bitscore:
                continue

            marker_id = parse_marker_id(sseqid)
            if marker_id not in marker_ids:
                continue

            marker_to_reads[marker_id].add(qseqid)
            cur = {"query": qseqid, "target": sseqid, "identity": pident, "aln_aa": alnlen, "bitscore": bitscore}
            prev = marker_best.get(marker_id)
            if prev is None or (cur["bitscore"], cur["identity"], cur["aln_aa"]) > (prev["bitscore"], prev["identity"], prev["aln_aa"]):
                marker_best[marker_id] = cur

    detail_rows = []
    detected = low_support = missing = 0
    pass_count = 0

    for marker_id in marker_ids:
        nreads = len(marker_to_reads.get(marker_id, set()))
        best = marker_best.get(marker_id)

        if best is None:
            missing += 1
            detail_rows.append({
                "marker_id": marker_id,
                "detected": 0,
                "unique_reads": 0,
                "best_identity": 0.0,
                "best_aln_aa": 0,
                "best_bitscore": 0.0,
                "best_query": "",
                "best_target": "",
                "status": "MISSING",
            })
            continue

        if nreads >= min_reads_per_marker:
            detected += 1
            pass_count += 1
            status = "PASS"
            detected_flag = 1
        else:
            low_support += 1
            status = "LOW_SUPPORT"
            detected_flag = 0

        detail_rows.append({
            "marker_id": marker_id,
            "detected": detected_flag,
            "unique_reads": nreads,
            "best_identity": round(best["identity"], 4),
            "best_aln_aa": best["aln_aa"],
            "best_bitscore": best["bitscore"],
            "best_query": best["query"],
            "best_target": best["target"],
            "status": status,
        })

    total_markers = len(marker_ids)
    c_pct = round(100.0 * detected / total_markers, 1) if total_markers else 0.0
    f_pct = round(100.0 * low_support / total_markers, 1) if total_markers else 0.0
    m_pct = round(100.0 * missing / total_markers, 1) if total_markers else 0.0
    gatekeeper_summary = f"FASTQ marker summary: C:{c_pct:.1f}%, F:{f_pct:.1f}%, M:{m_pct:.1f}%, n:{total_markers}"
    qc_score = round(100.0 * (pass_count / total_markers), 2) if total_markers else 0.0

    with open(detail_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=[
            "marker_id", "detected", "unique_reads",
            "best_identity", "best_aln_aa", "best_bitscore", "best_query", "best_target",
            "status"
        ])
        w.writeheader()
        w.writerows(detail_rows)

    summary = {
        "sample": args.sample,
        "platform": args.platform,
        "total_reads": total_reads,
        "detected_markers": detected,
        "total_markers": total_markers,
        "missing_markers": missing,
        "low_support_markers": low_support,
        "core_fraction": round(detected / total_markers, 4) if total_markers else 0.0,
        "required_markers": args.required_genes,
        "min_pident": min_pident,
        "min_aln_aa": min_aln_aa,
        "min_bitscore": min_bitscore,
        "min_reads_per_marker": min_reads_per_marker,
        "low_data_reads": low_data_reads,
        "max_target_seqs": max_target_seqs,
        "target_reporting_policy": "ALL_TARGETS" if max_target_seqs == 0 else "TOP_N_TARGETS",
        "long_read_target_censor_guard": ("PASS_ALL_TARGETS" if max_target_seqs == 0 else "FAIL_FINITE_TARGET_CAP") if long_read_profile else "NOT_APPLICABLE_SHORT_READ",
        "reported_hsps": reported_hsps,
        "queries_with_reported_hits": len(query_report_counts),
        "queries_at_report_limit": sum(1 for n in query_report_counts.values() if max_target_seqs > 0 and n >= max_target_seqs),
        "query_limit_saturation_fraction": round(sum(1 for n in query_report_counts.values() if max_target_seqs > 0 and n >= max_target_seqs) / len(query_report_counts), 4) if query_report_counts and max_target_seqs > 0 else 0.0,
        "threshold_status": defaults["threshold_status"],
        "engine_revision": ENGINE_REVISION,
        "gatekeeper_summary": gatekeeper_summary,
        "qc_score": qc_score,
        "pass_count": pass_count,
        "missing_count": missing,
        "low_support_count": low_support,
        "QC_STATUS": decide_marker_qc_status(
            detected, args.required_genes
        ),
        "LEGACY_QC_STATUS": decide_legacy_qc_status(
            detected, args.required_genes, total_reads, low_data_reads
        ),
        "DATA_SUFFICIENCY_STATUS": decide_data_sufficiency_status(
            total_reads, low_data_reads
        ),
        "data_sufficiency_basis": (
            "LEGACY_READ_COUNT_THRESHOLD_AUDIT_ONLY_NOT_QC_OVERRIDE"
        ),
    }

    with open(summary_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(summary.keys()))
        w.writeheader()
        w.writerow(summary)

    print("[DONE]")
    print("summary:", summary_csv)
    print("detail :", detail_csv)
    print("m8     :", m8)

if __name__ == "__main__":
    main()
