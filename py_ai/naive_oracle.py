"""
Brute-force Voronoi oracle. Assigns every pixel of a width x height grid to
its nearest site under the supplied distance callback. This is
O(width * height * sites) and is kept only as a test reference; the ported
library uses the Lee & Wong divide-and-conquer generator instead.
"""

import math
import random


def generate_voronoi_points(points, width, height, distance_callback):
    """
    Generate Voronoi points via a basic, naive algorithm. Takes any distance callback

    :param points: list of points [[x,y], ...]
    :param width: integer width
    :param height: integer height
    :param distance_callback: function(point1, point2) returning distance
    :returns: list of color arrays
    """
    colors = [
        {"point": e, "color": [math.ceil(random.random() * 255) for _ in range(3)]}
        for e in points
    ]

    image_data = []
    for index in range(width * height):
        coordinate = [index % height, math.ceil(index / height)]
        closest = colors[0]

        for c in colors:
            if isinstance(closest, list):
                # This branch unlikely in Python port
                pass
            elif distance_callback(closest["point"], coordinate) == distance_callback(
                c["point"], coordinate
            ):
                closest = [closest, c]
            else:
                closest = (
                    c
                    if distance_callback(c["point"], coordinate)
                    < distance_callback(closest["point"], coordinate)
                    else closest
                )

        # Find actual closest point
        best = colors[0]
        for c in colors:
            if distance_callback(c["point"], coordinate) < distance_callback(
                best["point"], coordinate
            ):
                best = c

        image_data.append(
            best["color"] if isinstance(best, dict) and "color" in best else [0, 0, 0]
        )

    return image_data
