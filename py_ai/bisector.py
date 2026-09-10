"""
Site and bisector graph helpers: object factories plus the predicates and
mutation helpers that keep the site <-> bisector graph consistent.
"""

from geometry import segment_intersection


def create_site(point):
    """Create a site object holding a raw coordinate."""
    return {"site": point, "bisectors": []}


def create_bisector(sites, up):
    """Create an empty bisector between two sites."""
    return {"sites": sites, "up": up, "points": [], "intersections": [], "compound": False}


def link_bisector_to_sites(bisector):
    """Register a bisector with both of its sites."""
    for site in bisector["sites"]:
        site["bisectors"].append(bisector)
    return bisector


def remove_bisector(bisector):
    """Unregister a bisector from both of its sites."""
    for site in bisector["sites"]:
        site["bisectors"] = [e for e in site["bisectors"] if e != bisector]


def find_hop_to(bisector, hop_from):
    """Find the other point across a bisector"""
    for e in bisector["sites"]:
        if e != hop_from:
            return e
    return bisector["sites"][0]


def is_point_on_edge(point, width, height):
    """Check if point is on an edge"""
    return point[0] == 0 or point[0] == width or point[1] == 0 or point[1] == height


def are_points_on_same_edge(p1, p2, width, height):
    """Check if two points are on the same edge"""
    return (
        (p1[0] == p2[0] and p1[0] == 0)
        or (p1[0] == p2[0] and p1[0] == width)
        or (p1[1] == p2[1] and p1[1] == 0)
        or (p1[1] == p2[1] and p1[1] == height)
    )


def is_bisector_trapped(trap_point, bisector, metric):
    """Determine if bisector is trapped in a site's polygon"""
    site0 = bisector["sites"][0]
    site1 = bisector["sites"][1]

    return all(
        metric["distance"](trap_point["site"], point) <= metric["distance"](site0["site"], point)
        and metric["distance"](trap_point["site"], point) <= metric["distance"](site1["site"], point)
        for point in bisector["points"]
    )


def get_extreme_point(bisector, go_up):
    """Find highest or lowest point of bisector"""
    if go_up:
        return max(point[1] for point in bisector["points"])
    else:
        return min(point[1] for point in bisector["points"])


def trim_bisector(target, intersector, intersection, metric):
    """Trim a bisector at a particular point"""
    if not intersector or not intersection:
        return

    # Find polygon site
    polygon_site = None
    for e in intersector["sites"]:
        if not any(d == e for d in target["sites"]):
            polygon_site = e
            break

    if not polygon_site:
        return

    # Filter points
    new_points = [
        point
        for point in target["points"]
        if metric["distance"](point, target["sites"][0]["site"])
        < metric["distance"](point, polygon_site["site"])
        and metric["distance"](point, target["sites"][1]["site"])
        < metric["distance"](point, polygon_site["site"])
    ]
    new_points.append(intersection)

    # Sort
    if target["up"]:
        target["points"] = sorted(new_points, key=lambda p: p[1])
    else:
        target["points"] = sorted(new_points, key=lambda p: p[0])


def bisector_intersection(b1, b2):
    """Find intersection of two bisectors"""
    if b1 == b2:
        return False

    for i in range(len(b1["points"]) - 1):
        for j in range(len(b2["points"]) - 1):
            intersect = segment_intersection(
                [b1["points"][i], b1["points"][i + 1]],
                [b2["points"][j], b2["points"][j + 1]],
            )
            if intersect:
                return intersect

    return False


def clear_out_orphans(orphanage, trap_point, metric):
    """Clear out orphans when a new merge line is created"""
    return [b for b in orphanage["bisectors"] if not is_bisector_trapped(trap_point, b, metric)]
