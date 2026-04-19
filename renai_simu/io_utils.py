"""I/O utilities: saving parameters, event logs, and summary statistics."""
from pathlib import Path
from datetime import datetime
from collections import Counter
import json
import numpy as np


def save_params_txt(sim, path):
    """Write human-readable parameters + final summary to `path`."""
    p = sim.p
    lines = []
    lines.append("=" * 64)
    lines.append("AffinityMD — Parameters & Final Summary")
    lines.append("=" * 64)
    lines.append(f"Run timestamp : {datetime.now().isoformat()}")
    lines.append("")

    lines.append("--- Parameters ---")
    for k in sorted(p.keys()):
        lines.append(f"  {k:20s} = {p[k]}")

    lines.append("")
    lines.append("--- Particles ---")
    for i in range(sim.N):
        typ = 'M' if sim.types[i] == 0 else 'F'
        lines.append(f"  {sim.labels[i]:14s} type={typ} "
                     f"age={sim.ages[i]:5.2f}")

    # Final state
    n_c   = int(np.sum(sim.couple[sim.males] != -1))
    max_c = min(sim.N_m, sim.N_f)

    lines.append("")
    lines.append("--- Final state ---")
    lines.append(f"Couples formed        : {n_c} / {max_c}")
    lines.append(f"Total stochastic events: {len(sim.event_log)}")

    lines.append("")
    lines.append("--- Final couples ---")
    if n_c == 0:
        lines.append("  (none)")
    for i in sim.males:
        j = int(sim.couple[i])
        if j != -1:
            lines.append(f"  {sim.labels[i]} <-> {sim.labels[j]}   "
                         f"B = {sim.B[i, j]:.3f}")

    lines.append("")
    lines.append("--- Singles ---")
    for i in range(sim.N):
        if sim.couple[i] == -1:
            typ = 'M' if sim.types[i] == 0 else 'F'
            lines.append(f"  {sim.labels[i]:14s} type={typ}")

    # Event breakdown
    ev_counts = Counter(k for (_, _, _, k) in sim.event_log)
    lines.append("")
    lines.append("--- Event breakdown ---")
    if not ev_counts:
        lines.append("  (no events)")
    for k, v in sorted(ev_counts.items()):
        lines.append(f"  {k:20s} : {v}")

    # Final B stats over M-F pairs only
    B_mf = sim.B[np.ix_(sim.males, sim.females)]
    lines.append("")
    lines.append("--- Final B_ij statistics (M-F only) ---")
    lines.append(f"  mean = {B_mf.mean():.3f}")
    lines.append(f"  std  = {B_mf.std():.3f}")
    lines.append(f"  min  = {B_mf.min():.3f}")
    lines.append(f"  max  = {B_mf.max():.3f}")
    lines.append(f"  # pairs with B > 0.5            : "
                 f"{int(np.sum(B_mf > 0.5))}")
    lines.append(f"  # pairs with B > B_thresh ({p['B_thresh']}) : "
                 f"{int(np.sum(B_mf > p['B_thresh']))}")

    Path(path).write_text('\n'.join(lines), encoding='utf-8')


def save_event_log(sim, path):
    """Save every stochastic event as TSV."""
    lines = ['time\tmale_idx\tfemale_idx\tmale_label\tfemale_label\tkind']
    for (t, i, j, kind) in sim.event_log:
        lines.append(
            f"{t:.3f}\t{i}\t{j}\t"
            f"{sim.labels[i]}\t{sim.labels[j]}\t{kind}"
        )
    Path(path).write_text('\n'.join(lines), encoding='utf-8')


def save_summary_json(sim, path):
    """Machine-readable JSON summary."""
    ev_counts = Counter(k for (_, _, _, k) in sim.event_log)
    B_mf = sim.B[np.ix_(sim.males, sim.females)]

    out = {
        'timestamp': datetime.now().isoformat(),
        'params':    sim.p,
        'particles': [
            {'index': int(i), 'label': sim.labels[i],
             'type':  'M' if sim.types[i] == 0 else 'F',
             'age':   float(sim.ages[i])}
            for i in range(sim.N)
        ],
        'final_couples': [
            {'male':   sim.labels[int(i)],
             'female': sim.labels[int(sim.couple[i])],
             'B':      float(sim.B[int(i), int(sim.couple[i])])}
            for i in sim.males if sim.couple[i] != -1
        ],
        'total_events':   len(sim.event_log),
        'event_breakdown': dict(ev_counts),
        'final_B_mf_stats': {
            'mean': float(B_mf.mean()),
            'std':  float(B_mf.std()),
            'min':  float(B_mf.min()),
            'max':  float(B_mf.max()),
        },
    }
    Path(path).write_text(json.dumps(out, indent=2), encoding='utf-8')


def save_all(sim, out_dir):
    """Write every output text/JSON file into `out_dir`."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    save_params_txt(sim,     out_dir / 'parmsdict.txt')
    save_event_log(sim,      out_dir / 'event_log.tsv')
    save_summary_json(sim,   out_dir / 'summary.json')
