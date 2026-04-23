import sys
import os
import random
import argparse

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from voronoi import generate_l1_voronoi


def get_color(n):
    colors = {
        4: "#4286f4",
        8: "#44f453",
        16: "#931d78",
        32: "#ff3c35",
        64: "#f4ad42",
        128: "#009182",
        256: "#993300",
        512: "#669999",
        1024: "#800000",
        2048: "#333300",
    }
    return colors.get(n, "#000000")


def hex_to_rgb(hex_color):
    hex_color = hex_color.lstrip("#")
    return tuple(int(hex_color[i : i + 2], 16) for i in (0, 2, 4))


def visualize(sites, width, height, show_bisectors=True, output_file=None):
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon as MplPolygon
    from matplotlib.collections import PatchCollection
    import numpy as np

    fig, ax = plt.subplots(figsize=(10, 10))
    ax.set_xlim(0, width)
    ax.set_ylim(0, height)
    ax.set_aspect("equal")
    ax.invert_yaxis()

    if show_bisectors:
        for site in sites:
            for bisector in site.get("bisectors", []):
                if bisector.get("points"):
                    xs = [p[0] for p in bisector["points"]]
                    ys = [p[1] for p in bisector["points"]]
                    merge_line = bisector.get("mergeLine")
                    color = get_color(merge_line) if merge_line else "#000000"
                    ax.plot(xs, ys, color=color, linewidth=0.5, alpha=0.6)

    patches = []
    colors = []
    for i, site in enumerate(sites):
        polygon_points = site.get("polygonPoints")
        if not polygon_points:
            continue

        polygon_points = polygon_points + [polygon_points[0]]
        poly = MplPolygon(polygon_points, closed=True)
        patches.append(poly)
        colors.append(i / len(sites))

    collection = PatchCollection(patches, cmap="tab20", alpha=0.3)
    collection.set_array(np.array(colors))
    ax.add_collection(collection)

    for site in sites:
        ax.plot(site["site"][0], site["site"][1], "k.", markersize=2)

    title = f"L1 Voronoi Diagram - {len(sites)} sites"
    if show_bisectors:
        title += " (with bisectors)"
    ax.set_title(title)
    ax.set_xlabel("x")
    ax.set_ylabel("y")

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=150, bbox_inches="tight")
        print(f"Saved to {output_file}")
    else:
        plt.show()

    plt.close()


def main():
    parser = argparse.ArgumentParser(description="L1 Voronoi Diagram Visualizer")
    parser.add_argument("-n", "--points", type=int, default=32, help="number of points")
    parser.add_argument("-w", "--width", type=int, default=400, help="canvas width")
    parser.add_argument("-H", "--height", type=int, default=400, help="canvas height")
    parser.add_argument("--no-bisectors", action="store_true", help="hide bisectors")
    parser.add_argument("-o", "--output", type=str, help="output file (PNG/PDF)")

    args = parser.parse_args()

    width = args.width
    height = args.height

    points = [
        [int(random.random() * width), int(random.random() * height)]
        for _ in range(args.points)
    ]

    print(f"Generating L1 Voronoi for {args.points} points...")
    print(f"Width: {width}, Height: {height}")

    sites = generate_l1_voronoi(points, width, height, True)

    print(f"Generated {len(sites)} site polygons")

    visualize(
        sites,
        width,
        height,
        show_bisectors=not args.no_bisectors,
        output_file=args.output,
    )


if __name__ == "__main__":
    main()
