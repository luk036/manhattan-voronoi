"""Python side of the JS <-> Python differential harness.

Reads the shared corpus JSON (written by tests/differential.js), runs each case
through the py_ai port with nudge disabled, and writes structural summaries that
differential.js compares against the JS summaries.
"""

import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "py_ai"))

from voronoi import generate_l1_voronoi  # noqa: E402


def summarize(c):
    sites = generate_l1_voronoi(
        [list(p) for p in c["sites"]], c["w"], c["h"], False
    )
    out = []
    for site in sites:
        neighbors = sorted(
            [list(n) for n in site.get("neighbors", [])], key=lambda p: (p[0], p[1])
        )
        polygon = site.get("polygon_points") if "polygon_points" in site else None
        out.append(
            {
                "site": list(site["site"]),
                "neighbors": neighbors,
                "polyCount": len(polygon) if polygon is not None else None,
            }
        )
    return out


def main():
    corpus_path, out_path = sys.argv[1], sys.argv[2]
    with open(corpus_path, "r", encoding="utf-8") as f:
        cases = json.load(f)
    results = [summarize(c) for c in cases]
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(results, f)


if __name__ == "__main__":
    main()
