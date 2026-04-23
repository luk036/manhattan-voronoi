import sys
import os

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


def random_normal(sharpness):
    return sum([__import__("random").random() for _ in range(sharpness)]) / sharpness


def main():
    width = 400
    height = 400

    num_points = 32

    import random

    raw = [
        [int(random.random() * width), int(random.random() * height)]
        for _ in range(num_points)
    ]
    sites = raw[:]

    print(f"Generating L1 Voronoi for {num_points} points...")
    print(f"Width: {width}, Height: {height}")

    vector_points = generate_l1_voronoi(sites, width, height, True)

    print(f"Generated {len(vector_points)} site polygons")

    for i, site in enumerate(vector_points[:3]):
        print(f"  Site {i}: site={site['site']}, bisectors={len(site['bisectors'])}")

    print("\nDemo completed successfully!")


if __name__ == "__main__":
    main()
