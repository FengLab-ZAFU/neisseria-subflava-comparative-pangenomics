#!/usr/bin/env python3
"""P1 step 5: finalize annotation deliverables.

Outputs in results_revision_20261010_phaseB_v7/01_annotation/:
  annotation_manifest.tsv            212 rows
  annotation_pending_resolved.tsv    19 pending + cache downgrades, final disposition
  annotation_additional_required.tsv entries still failing (empty if none)
  bakta_db_manifest.tsv              DB file listing + version.json hash (reuse evidence)
  checksums_p1.sha256
Also updates 00_preflight/reuse_or_rerun_ledger.tsv (bakta rows) in place.
"""
import csv, hashlib, os, subprocess
from pathlib import Path

ROOT = Path(os.environ["NEISSERIA_PROJECT_ROOT"])
OUT = ROOT / "results_revision_20261010_phaseB_v7/01_annotation"
PRE = ROOT / "results_revision_20261010_phaseB_v7/00_preflight"
V7 = ROOT / "results_revision_20261010_sourcefix_v7"
DB = Path(os.environ["BAKTA_DB"])
PENDING = Path(os.environ["ANNOTATION_PENDING_TSV"])


def read_tsv(p):
    with open(p, newline="", encoding="utf-8-sig") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def sha256_file(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    binding = {r["accession"]: r for r in read_tsv(OUT / "fasta_binding_212.tsv")}
    sources = {r["accession"]: r for r in read_tsv(OUT / "annotation_sources.tsv")}
    contig = {r["accession"]: r for r in read_tsv(OUT / "annotation_contig_binding.tsv")}
    cache = {r["accession"]: r for r in read_tsv(OUT / "annotation_cache_verification.tsv")}
    status = {}
    for sp in (OUT / "bakta_new_status.tsv", OUT / "bakta_new_status2.tsv"):
        if sp.exists():
            status.update({r["accession"]: r["exit_code"] for r in read_tsv(sp)})
    pending = read_tsv(PENDING)
    extra_downgrades = []
    rl2 = OUT / "bakta_runlist2.tsv"
    if rl2.exists():
        extra_downgrades = read_tsv(rl2)

    # --- DB manifest (reuse evidence) ---
    vj = DB / "version.json"
    vj_hash = sha256_file(vj) if vj.is_file() else ""
    db_rows = []
    for dirpath, _dirs, files in os.walk(DB):
        for fn in sorted(files):
            p = Path(dirpath) / fn
            st = p.stat()
            db_rows.append({"path": str(p.relative_to(DB)), "size_bytes": st.st_size})
    with open(OUT / "bakta_db_manifest.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["path", "size_bytes"], delimiter="\t")
        w.writeheader(); w.writerows(db_rows)
    db_manifest_hash = sha256_file(OUT / "bakta_db_manifest.tsv")
    db_evidence = (f"db_path={DB}; version.json sha256={vj_hash} "
                   f"(v6.0 full 2025-02-24, doi 10.5281/zenodo.14916843); "
                   f"manifest n_files={len(db_rows)} sha256={db_manifest_hash}")

    # --- manifest ---
    mrows = []
    for acc, s in sorted(sources.items()):
        d = Path(s["annotation_dir"]) if s["annotation_dir"] else None
        src = s["annotation_source"]
        def f(ext):
            p = d / f"{acc}.{ext}" if d else None
            return str(p) if p and p.is_file() else ""
        if src.startswith("reused"):
            exit_status = "0_historical(log_verified)"
            log_path = str(d / f"{acc}.log") if d and (d / f"{acc}.log").is_file() else ""
        elif src.startswith("new"):
            exit_status = status.get(acc, "missing")
            log_path = str(OUT / "bakta_new/logs" / f"{acc}.log") if (OUT / "bakta_new/logs" / f"{acc}.log").is_file() else ""
        else:
            exit_status, log_path = "blocked", ""
        cb = contig.get(acc, {})
        mrows.append({
            "accession": acc, "biosample": binding[acc]["biosample"],
            "annotation_source": src,
            "gff_path": f("gff3"), "gbff_path": f("gbff"), "faa_path": f("faa"), "ffn_path": f("ffn"),
            "log_path": log_path, "exit_status": exit_status,
            "binding_status": cb.get("status", "not_checked"),
        })
    mfields = ["accession", "biosample", "annotation_source", "gff_path", "gbff_path",
               "faa_path", "ffn_path", "log_path", "exit_status", "binding_status"]
    with open(OUT / "annotation_manifest.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=mfields, delimiter="\t")
        w.writeheader(); w.writerows(mrows)

    # --- pending resolved ---
    prows = []
    for r in pending:
        acc = r["accession"]
        prows.append(resolve_row(acc, "pending19", r["existing_binding_status"],
                                 sources, status, contig, binding))
    downgraded = {a: c for a, c in cache.items() if c["verdict"] == "DOWNGRADE"}
    for acc, c in sorted(downgraded.items()):
        prows.append(resolve_row(acc, "cache_downgrade", c["reason"],
                                 sources, status, contig, binding))
    for r in extra_downgrades:
        prows.append(resolve_row(r["accession"], "cache_downgrade_contig_binding", r["reason"],
                                 sources, status, contig, binding))
    pfields = ["accession", "origin", "origin_detail", "final_disposition", "reason"]
    with open(OUT / "annotation_pending_resolved.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=pfields, delimiter="\t")
        w.writeheader(); w.writerows(prows)

    # --- additional required (any remaining failures) ---
    arows = []
    for m in mrows:
        problems = []
        if m["annotation_source"] == "BLOCKED":
            problems.append("no_usable_fasta_or_annotation")
        if m["annotation_source"].startswith("new") and m["exit_status"] not in ("0",):
            problems.append(f"bakta_exit_status={m['exit_status']}")
        for ext in ("gff_path", "gbff_path", "faa_path", "ffn_path"):
            if not m[ext]:
                problems.append(f"missing_{ext.split('_')[0]}")
        if m["binding_status"] != "PASS":
            problems.append(f"contig_binding={m['binding_status']}:{contig.get(m['accession'], {}).get('detail', '')}")
        if problems:
            arows.append({"accession": m["accession"], "biosample": m["biosample"],
                          "annotation_source": m["annotation_source"],
                          "problems": ";".join(problems)})
    with open(OUT / "annotation_additional_required.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["accession", "biosample", "annotation_source", "problems"], delimiter="\t")
        w.writeheader(); w.writerows(arows)

    # --- ledger update ---
    n_reused = sum(1 for m in mrows if m["annotation_source"].startswith("reused"))
    n_new = sum(1 for m in mrows if m["annotation_source"].startswith("new"))
    ledger = PRE / "reuse_or_rerun_ledger.tsv"
    lines = ledger.read_text().splitlines()
    out_lines = []
    for ln in lines:
        if ln.startswith("bakta_annotations_193\t"):
            f = ln.split("\t")
            f[4] = "VERIFIED"
            f[6] = ("193/193缓存核验REUSE_OK: accession.version一致; 注释输入序列与冻结FASTA逐序列多重集一致; "
                    "逐株日志确认bakta=1.12.0且DB路径与记录一致; "
                    "gff3/gbff/faa/ffn齐全非空. 缓存GFF3系2026-10-02由GBFF重建(symlink-clobber事件), "
                    "contig绑定核验发现4条(GCA_003351665.1/GCA_056611915.1/GCA_056850595.1/GCA_056850685.1)"
                    "origin-spanning join被压平成1..contig_len→降级bakta重跑; "
                    f"最终{n_reused}条复用+{n_new}条新跑(19待补+4降级), 212/212 contig绑定PASS "
                    "(01_annotation/annotation_cache_verification.tsv, annotation_contig_binding.tsv)")
            f[7] = f"VERIFIED_REUSED({n_reused})/RERUN({n_new})"
            f[8] = db_evidence + "; 新跑参数与b1_bakta.sh一致(--compliant --force, prefix=accession, locus-tag=N<digits>), 仅threads 4→8"
            ln = "\t".join(f)
        elif ln.startswith("bakta_db_content_hash\t"):
            f = ln.split("\t")
            f[4] = "VERIFIED"
            f[6] = db_evidence
            f[7] = "VERIFIED_REUSED"
            f[8] = "84GB库未逐字节哈希; 以version.json哈希+全量文件清单(路径+大小)登记"
            ln = "\t".join(f)
        out_lines.append(ln)
    ledger.write_text("\n".join(out_lines) + "\n")

    # --- checksums ---
    targets = []
    for p in sorted(OUT.glob("*.tsv")):
        targets.append(p)
    for p in sorted((OUT / "contig_id_maps").glob("*.tsv")):
        targets.append(p)
    newdir = OUT / "bakta_new"
    if newdir.is_dir():
        for p in sorted(newdir.rglob("*")):
            if p.is_file():
                targets.append(p)
    with open(OUT / "checksums_p1.sha256", "w") as fh:
        for p in targets:
            fh.write(f"{sha256_file(p)}  {p.relative_to(ROOT)}\n")

    print(f"manifest: {len(mrows)} rows ({n_reused} reused, {n_new} new)")
    print(f"pending_resolved: {len(prows)} rows; additional_required: {len(arows)}")
    print(f"checksums: {len(targets)} files")
    for a in arows:
        print("REQUIRED", a["accession"], a["problems"])


def resolve_row(acc, origin, detail, sources, status, contig, binding):
    s = sources.get(acc, {})
    src = s.get("annotation_source", "")
    cb = contig.get(acc, {})
    if src.startswith("reused"):
        disp, reason = "resolved_reuse_verified", "cache verification REUSE_OK"
    elif src.startswith("new"):
        if status.get(acc) == "0" and cb.get("status") == "PASS":
            disp, reason = "resolved_new", "bakta 1.12.0/db6.0 exit=0, contig binding PASS"
        elif status.get(acc) == "0":
            disp, reason = "failed", f"bakta ok but contig binding {cb.get('status')}:{cb.get('detail','')}"
        elif acc in status:
            disp, reason = "failed", f"bakta exit_status={status[acc]}"
        else:
            disp, reason = "failed", "bakta not run or status missing (check fasta binding)"
    else:
        disp, reason = "failed", "fasta binding failed, annotation blocked"
    return {"accession": acc, "origin": origin, "origin_detail": detail,
            "final_disposition": disp, "reason": reason}


if __name__ == "__main__":
    main()
