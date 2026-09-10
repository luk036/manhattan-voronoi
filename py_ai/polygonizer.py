"""
Post-processing stage: for each site, chain its bisectors into an ordered
polygon, then enrich the site with the polygon points and its neighbor list.
"""

from geometry import angle, same_point
from bisector import (
    is_point_on_edge,
    are_points_on_same_edge,
    bisector_intersection,
    find_hop_to,
)


def chain_bisector_points(site, metric):
    """Chain the bisectors of a site into a closed polygon point list."""
    bisectors = site["bisectors"]

    if not bisectors:
        return None

    # Find starting bisector
    start_bisector = None
    for b in bisectors:
        if any(is_point_on_edge(p, metric["width"], metric["height"]) for p in b["points"]):
            start_bisector = b
            break
    start_bisector = start_bisector or bisectors[0]

    starting_points = start_bisector["points"]

    # Reverse if ends on edge
    if is_point_on_edge(starting_points[-1], metric["width"], metric["height"]):
        starting_points = starting_points[::-1]

    polygon_points = list(starting_points)
    used = [start_bisector]

    # Walk through remaining bisectors
    for _ in range(len(bisectors) - 1):
        last = polygon_points[-1]

        # Find next bisector with closest endpoint
        next_bisector = None
        min_dist = float("inf")

        for b in bisectors:
            if b in used:
                continue

            d1 = metric["distance"](last, b["points"][0])
            d2 = metric["distance"](last, b["points"][-1])
            d = min(d1, d2)

            if d < min_dist:
                min_dist = d
                next_bisector = b

        if next_bisector:
            next_points = next_bisector["points"]
            if same_point(next_points[-1], last):
                next_points = next_points[::-1]

            polygon_points.extend(next_points)
            used.append(next_bisector)

    return polygon_points


def append_open_edge_corners(site, polygon_points, metric):
    """Fill in canvas corners when the polygon opens onto two different edges."""
    if polygon_points is None:
        return polygon_points

    corners = [
        [0, 0],
        [metric["width"], 0],
        [metric["width"], metric["height"]],
        [0, metric["height"]],
    ]

    # Handle case where polygon starts and ends on different edges
    if is_point_on_edge(
        polygon_points[0], metric["width"], metric["height"]
    ) and is_point_on_edge(polygon_points[-1], metric["width"], metric["height"]):
        if not are_points_on_same_edge(
            polygon_points[0], polygon_points[-1], metric["width"], metric["height"]
        ):
            filtered_corners = [
                e
                for e in corners
                if not any(
                    bisector_intersection({"points": [e, site["site"]]}, d)
                    for d in site["bisectors"]
                )
            ]
            polygon_points.extend(filtered_corners)

    return polygon_points


def polygonize_site(site, metric):
    """Polygonize and enrich a site in place: adds polygon_points and neighbors."""
    polygon_points = chain_bisector_points(site, metric)

    if polygon_points is None:
        return site

    polygon_points = append_open_edge_corners(site, polygon_points, metric)

    # Sort by angle
    site["polygon_points"] = sorted(
        polygon_points, key=lambda p: angle(site["site"], p)
    )

    # Get neighbors
    site["neighbors"] = [find_hop_to(b, site)["site"] for b in site["bisectors"]]

    return site
