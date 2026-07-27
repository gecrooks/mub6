"""Merge local + Modal X-cover results, tally coverage, list gap lines
for the continuation patch pass."""

import re
import sys

def main(modal_log="run_modal_cover.log", local_log="xcover_lines.log"):
    lines = {}
    try:
        for ln in open(local_log):
            m = re.match(r"c=([+-][\d.]+) t=\[([+-][\d.]+),([+-][\d.]+)\] "
                         r"anchors=(\d+) gaps=(\d+)", ln)
            if m:
                lines[round(float(m.group(1)), 4)] = dict(
                    t_lo=float(m.group(2)), t_hi=float(m.group(3)),
                    anchors=int(m.group(4)), gaps=int(m.group(5)), src="local")
    except FileNotFoundError:
        pass
    pat = re.compile(r"'c': ([+-]?[\d.]+).*?'t_lo': np\.float64\(([-\d.e]+)\)"
                     r".*?'t_hi': np\.float64\(([-\d.e]+)\).*?"
                     r"'anchors': (\d+), 'gaps': (\d+)")
    for ln in open(modal_log):
        m = pat.search(ln)
        if m:
            lines[round(float(m.group(1)), 4)] = dict(
                t_lo=float(m.group(2)), t_hi=float(m.group(3)),
                anchors=int(m.group(4)), gaps=int(m.group(5)), src="modal")
    tot_a = sum(v["anchors"] for v in lines.values())
    tot_g = sum(v["gaps"] for v in lines.values())
    gap_lines = sorted(c for c, v in lines.items() if v["gaps"] > 0)
    width = sum(v["t_hi"] - v["t_lo"] for v in lines.values())
    print(f"lines: {len(lines)}  anchors: {tot_a}  gaps: {tot_g}")
    print(f"summed line-width: {width:.2f} (region area covered ~ "
          f"{width * 0.025:.3f})")
    print(f"gap lines ({len(gap_lines)}): {gap_lines}")
    n_modal = sum(1 for v in lines.values() if v["src"] == "modal")
    print(f"sources: {n_modal} modal, {len(lines) - n_modal} local")


if __name__ == "__main__":
    main(*sys.argv[1:])
