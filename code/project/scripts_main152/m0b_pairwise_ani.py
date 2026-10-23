#!/usr/bin/env python3
"""M0b: verify the 152x152 pairwise FastANI cache against frozen FASTA hashes
and emit directional/unordered tables + symmetric matrices.

Rules (plan M0b + config):
- reuse cache only if every cached input FASTA sha256 == frozen sha256 (152/152)
- 23,104 directional pairs incl. 152 self pairs; keep measured self ANI
- unordered ANI = mean(ani_ab, ani_ba); AF_min = min(af_ab, af_ba)
- 0<=ANI<=100, 0<=AF<=1; no unexplained NA; matrix symmetric by construction
- no re-screening by ANI>=95
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import RESULTS, PACKET, read_roster, read_tsv, write_tsv

OUT = os.path.join(RESULTS, "00_preflight/pairwise_ANI")
REF = os.path.join(PACKET, "00_input/reference_sources")
IN_CHK = os.path.join(RESULTS, "00_preflight/input_check")


def main():
    roster = read_roster()
    binding = {r["accession"]: r for r in read_tsv(os.path.join(IN_CHK, "input_binding.tsv"))}

    # 1) cache input binding vs actual frozen FASTA hashes
    cache_bind = {r["accession"]: r for r in read_tsv(os.path.join(REF, "pairwise_ANI_cached_input_binding_main152.tsv"))}
    audit_rows, mism = [], []
    for acc in roster:
        cb = cache_bind.get(acc)
        if cb is None:
            mism.append((acc, "missing from cache binding")); continue
        ok = cb["sha256"] == binding[acc]["actual_sha256"]
        if not ok:
            mism.append((acc, f"cache {cb['sha256'][:12]} != actual {binding[acc]['actual_sha256'][:12]}"))
        audit_rows.append(dict(accession=acc, cached_sha256=cb["sha256"],
                               actual_frozen_sha256=binding[acc]["actual_sha256"],
                               match="yes" if ok else "NO",
                               cache_unit_id=cb["unit_id"]))
    write_tsv(os.path.join(OUT, "cache_reuse_audit.tsv"),
              ["accession", "cached_sha256", "actual_frozen_sha256", "match", "cache_unit_id"],
              audit_rows)
    if mism:
        print("CACHE HASH MISMATCH — rerun FastANI required for:", mism); sys.exit(2)
    print("cache binding: 152/152 sha256 match -> reuse authorized")

    # 2) directional pairs validation
    raw = read_tsv(os.path.join(REF, "pairwise_ANI_cached_raw_main152.tsv"))
    assert len(raw) == 23104, f"directional rows {len(raw)} != 23104"
    seen = {}
    bad = []
    for r in raw:
        q, ref = r["query_accession"], r["reference_accession"]
        ani = float(r["ANI_percent"]); mf = int(r["matched_fragments"]); tf = int(r["total_fragments"])
        if q not in binding or ref not in binding:
            bad.append((q, ref, "accession not in roster")); continue
        if not (0 <= ani <= 100) or mf > tf or tf <= 0:
            bad.append((q, ref, f"range ani={ani} mf={mf} tf={tf}")); continue
        seen[(q, ref)] = (ani, mf / tf)
    if len(seen) != 23104:
        dupes = 23104 - len(seen)
        bad.append(("GLOBAL", "GLOBAL", f"duplicate directional pairs: {dupes}"))
    if bad:
        print("BAD ROWS:", bad[:10]); sys.exit(3)

    write_tsv(os.path.join(OUT, "directional_pairs.tsv"),
              ["query_accession", "reference_accession", "ANI_percent", "AF"],
              [dict(query_accession=q, reference_accession=r, ANI_percent=f"{v[0]:.6f}".rstrip("0").rstrip("."),
                    AF=f"{v[1]:.6f}") for (q, r), v in sorted(seen.items())])

    # 3) unordered pairs + matrices
    un_rows = []
    ani_mat = {a: {} for a in roster}
    af_mat = {a: {} for a in roster}
    for i, a in enumerate(roster):
        for b in roster[i + 1:]:
            ani_ab, af_ab = seen[(a, b)]
            ani_ba, af_ba = seen[(b, a)]
            ani_mean = (ani_ab + ani_ba) / 2
            af_min = min(af_ab, af_ba)
            un_rows.append(dict(accession_a=a, accession_b=b,
                                ani_ab=f"{ani_ab:.6f}".rstrip("0").rstrip("."),
                                ani_ba=f"{ani_ba:.6f}".rstrip("0").rstrip("."),
                                ani_mean=f"{ani_mean:.6f}",
                                af_ab=f"{af_ab:.6f}", af_ba=f"{af_ba:.6f}", af_min=f"{af_min:.6f}"))
            ani_mat[a][b] = ani_mat[b][a] = ani_mean
            af_mat[a][b] = af_mat[b][a] = af_min
    assert len(un_rows) == 11476
    for a in roster:
        ani_mat[a][a] = seen[(a, a)][0]   # measured self ANI, not forced 100
        af_mat[a][a] = seen[(a, a)][1]
    write_tsv(os.path.join(OUT, "unordered_pairs.tsv"),
              ["accession_a", "accession_b", "ani_ab", "ani_ba", "ani_mean", "af_ab", "af_ba", "af_min"],
              un_rows)

    for name, mat, fmt in [("ani_mean_matrix.tsv", ani_mat, "{:.6f}"),
                           ("af_min_matrix.tsv", af_mat, "{:.6f}")]:
        with open(os.path.join(OUT, name), "w") as fh:
            fh.write("accession\t" + "\t".join(roster) + "\n")
            for a in roster:
                fh.write(a + "\t" + "\t".join(fmt.format(mat[a][b]) for b in roster) + "\n")

    # input binding copy for the ANI dataset
    write_tsv(os.path.join(OUT, "input_binding.tsv"),
              ["accession", "unit_id", "fasta_path", "sha256", "cache_reused"],
              [dict(accession=acc, unit_id=cache_bind[acc]["unit_id"],
                    fasta_path=binding[acc]["path"], sha256=binding[acc]["actual_sha256"],
                    cache_reused="yes_hash_verified") for acc in roster])

    self_anis = [seen[(a, a)][0] for a in roster]
    off = [r["ani_mean"] for r in un_rows]
    summary = dict(
        module="M0b", status="PASS",
        cache_binding="152/152 sha256 match",
        directional_pairs=len(seen), unordered_nonself_pairs=len(un_rows),
        self_ani_min=min(self_anis), self_ani_max=max(self_anis),
        offdiag_ani_min=min(float(x) for x in off), offdiag_ani_max=max(float(x) for x in off),
        matrix_order="main152_accessions.txt row order",
        note="symmetric matrices built from mean/min of both directions; self pairs keep measured values; "
             "no ANI>=95 rescreening applied",
    )
    with open(os.path.join(OUT, "pairwise_ani_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=1)
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
