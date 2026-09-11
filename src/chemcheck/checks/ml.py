"""ML-readiness checks: measurements, target shift, leakage suspects, outliers."""
from __future__ import annotations

from collections import defaultdict

from ..models import Finding, Severity

MAX_ROWS = 500


def _labels(records):
    return [(r, r.label) for r in records if r.valid and r.label is not None
            and str(r.label) not in ("", "nan", "None")]


def _is_numeric(vals) -> bool:
    try:
        [float(v) for v in vals]
        return True
    except Exception:
        return False


def check_duplicated_measurements(records, ctx):
    g = defaultdict(list)
    for r in records:
        if r.valid and r.canon_smi and r.label is not None and str(r.label) != "nan":
            g[(r.canon_smi, str(r.label))].append(r.row_id)
    dups = {k: v for k, v in g.items() if len(v) > 1}
    if not dups:
        return None
    affected = sorted({rid for v in dups.values() for rid in v})
    ex = [{"structure": k[0][:80], "label": k[1], "rows": v[:6]} for k, v in
          sorted(dups.items(), key=lambda kv: -len(kv[1]))[:ctx["max_examples"]]]
    return Finding(check_id="duplicated_measurements", severity=Severity.WARNING,
                   title="duplicated measurements (same structure, same label)",
                   count=len(affected), rate=len(affected) / max(1, len(records)),
                   affected_rows=affected[:MAX_ROWS], total_affected=len(affected),
                   examples=ex,
                   recommendation="Aggregate replicates (mean/median) instead of repeating rows — duplicates "
                                  "overweight those points and leak across splits.",
                   details=f"{len(dups)} (structure, label) groups repeated.")


def check_conflicting_measurements(records, ctx):
    g = defaultdict(list)
    for r in records:
        if r.valid and r.canon_smi and r.label is not None and str(r.label) not in ("", "nan", "None"):
            g[r.canon_smi].append((r.row_id, r.label))
    conflicts = {}
    for k, v in g.items():
        labs = [x[1] for x in v]
        if len(set(str(x) for x in labs)) <= 1:
            continue
        if _is_numeric(labs):
            f = [float(x) for x in labs]
            spread = max(f) - min(f)
            # conflict if range > 1.0 log unit (or >10% of range heuristic thin)
            if spread > 1.0:
                conflicts[k] = (v, spread)
        else:
            conflicts[k] = (v, None)
    if not conflicts:
        return None
    affected = sorted({rid for v, _ in conflicts.values() for rid, _ in v})
    ex = []
    for k, (v, spread) in sorted(conflicts.items(), key=lambda kv: -(kv[1][1] or 99))[:ctx["max_examples"]]:
        ex.append({"structure": k[:80], "rows": [rid for rid, _ in v][:6],
                   "labels": [str(l) for _, l in v][:6],
                   "spread": round(spread, 3) if spread is not None else "categorical-mismatch"})
    return Finding(check_id="conflicting_measurements", severity=Severity.ERROR,
                   title="conflicting measurements (same structure, different labels)",
                   count=len(affected), rate=len(affected) / max(1, len(records)),
                   affected_rows=affected[:MAX_ROWS], total_affected=len(affected),
                   examples=ex,
                   recommendation="Curate conflicts: check assay conditions, average with uncertainty, or drop "
                                  "the structure. Models cannot learn contradictory labels.",
                   details=f"{len(conflicts)} structures with incompatible labels.")


def check_target_shift(records, ctx):
    tr = [r for r in records if r.valid and (r.split or "").lower() == "train"]
    te = [r for r in records if r.valid and (r.split or "").lower() in ("test", "valid", "validation")]
    if not tr or not te:
        return None
    tr_labs = [r.label for r in tr if r.label is not None and str(r.label) not in ("", "nan", "None")]
    te_labs = [r.label for r in te if r.label is not None and str(r.label) not in ("", "nan", "None")]
    if len(tr_labs) < 10 or len(te_labs) < 5:
        return None
    if _is_numeric(tr_labs) and _is_numeric(te_labs):
        import numpy as np
        a = np.array([float(x) for x in tr_labs], dtype=float)
        b = np.array([float(x) for x in te_labs], dtype=float)
        try:
            from scipy.stats import ks_2samp
            stat, p = ks_2samp(a, b)
        except Exception:
            stat, p = 0.0, 1.0
        mean_shift = abs(float(a.mean()) - float(b.mean()))
        pooled = float(a.std() + b.std()) / 2 or 1.0
        d = mean_shift / pooled
        if p >= 0.05 and d < 0.5:
            return None
        sev = Severity.WARNING if d < 0.8 else Severity.ERROR
        return Finding(check_id="target_distribution_shift", severity=sev,
                       title=f"target distribution shift (train mean={a.mean():.2f} vs test mean={b.mean():.2f})",
                       count=len(te_labs), rate=0.0, affected_rows=[r.row_id for r in te][:MAX_ROWS],
                       total_affected=len(te_labs),
                       examples=[{"train_mean": round(float(a.mean()), 3),
                                  "test_mean": round(float(b.mean()), 3),
                                  "ks_stat": round(float(stat), 3), "ks_p": float(f"{p:.3g}"),
                                  "cohens_d": round(float(d), 3)}],
                       recommendation="Stratify splits by label and report per-split histograms; a shifted test "
                                      "set measures extrapolation, not interpolation — label it as such.",
                       details=f"KS p={p:.3g}, Cohen's d={d:.2f}.")
    # categorical
    from collections import Counter
    ca, cb = Counter(str(x) for x in tr_labs), Counter(str(x) for x in te_labs)
    Ma = sum(ca.values()) or 1
    Mb = sum(cb.values()) or 1
    classes = set(ca) | set(cb)
    diff = max(abs(ca[c] / Ma - cb[c] / Mb) for c in classes)
    if diff < 0.15:
        return None
    return Finding(check_id="target_distribution_shift", severity=Severity.WARNING,
                   title=f"label prevalence shift (max Δ={diff * 100:.1f}pp)",
                   count=len(te_labs), rate=0.0, affected_rows=[r.row_id for r in te][:MAX_ROWS],
                   total_affected=len(te_labs),
                   examples=[{"class": c, "train_frac": round(ca[c] / Ma, 3),
                              "test_frac": round(cb[c] / Mb, 3)} for c in sorted(classes)][:ctx["max_examples"]],
                   recommendation="Stratify by class; report prevalence per split.",
                   details="")


def check_outliers(records, ctx):
    lab = _labels(records)
    if len(lab) < 10:
        return None
    vals = [lbl for _, lbl in lab]
    if not _is_numeric(vals):
        return None
    import numpy as np
    a = np.array([float(v) for v in vals])
    mu, sd = float(a.mean()), float(a.std()) or 1.0
    out = [(r.row_id, float(r.label)) for r, _ in lab if abs(float(r.label) - mu) / sd > 4]
    if not out:
        return None
    ex = [{"row": rid, "label": v, "z": round((v - mu) / sd, 2)} for rid, v in out[:ctx["max_examples"]]]
    return Finding(check_id="label_outliers", severity=Severity.INFO,
                   title=f"{len(out)} label outliers (|z|>4)",
                   count=len(out), rate=len(out) / max(1, len(records)),
                   affected_rows=[rid for rid, _ in out][:MAX_ROWS], total_affected=len(out),
                   examples=ex,
                   recommendation="Verify units (nM vs µM is a classic) and assay ceilings; winsorize or model "
                                  "with robust losses rather than silently keeping them.",
                   details=f"mean={mu:.3f} sd={sd:.3f}.")


def check_split_label_leakage(records, ctx):
    # thin "target leakage": does split membership predict the label?
    tr = [r for r in records if r.valid and (r.split or "").lower() == "train"]
    te = [r for r in records if r.valid and (r.split or "").lower() in ("test", "valid", "validation")]
    if not tr or not te:
        return None
    tr_labs = [r.label for r in tr if r.label is not None and str(r.label) not in ("", "nan", "None")]
    te_labs = [r.label for r in te if r.label is not None and str(r.label) not in ("", "nan", "None")]
    if len(tr_labs) < 10 or len(te_labs) < 5 or not _is_numeric(tr_labs + te_labs):
        return None
    import numpy as np
    a = np.array([float(x) for x in tr_labs])
    b = np.array([float(x) for x in te_labs])
    d = abs(float(a.mean() - b.mean())) / ((float(a.std() + b.std())) / 2 or 1.0)
    if d < 0.8:
        return None
    return Finding(check_id="target_leakage_split_predicts_label", severity=Severity.WARNING,
                   title="split predicts target (possible target/provenance leakage)",
                   count=len(te_labs), rate=0.0, affected_rows=[r.row_id for r in te][:MAX_ROWS],
                   total_affected=len(te_labs),
                   examples=[{"train_mean": round(float(a.mean()), 3), "test_mean": round(float(b.mean()), 3),
                              "cohens_d": round(float(d), 3)}],
                   recommendation="Check for provenance leakage (e.g. one assay/series per split, time trends, "
                                  "or label-derived features). Re-split stratified by source and time.",
                   details="Large label gap between splits often means the split, not chemistry, explains performance.")


CHECKS = [check_duplicated_measurements, check_conflicting_measurements,
          check_target_shift, check_outliers, check_split_label_leakage]
