"""
The L1 metric: bisector construction and the upward-test used while walking
merge lines. The rest of the pipeline consumes these through the object
returned by `create_l1_metric`, so this module is the seam where a different
metric (e.g. L-infinity) would slot in.
"""

from geometry import same_point, distance
from bisector import create_bisector


def create_l1_metric(width, height):
    """Build the L1 metric strategy for a fixed canvas."""

    def bisector(p1, p2):
        return find_l1_bisector(p1, p2, width, height)

    return {
        "width": width,
        "height": height,
        "distance": distance,
        "bisector": bisector,
        "is_upward": is_new_bisector_upward,
    }


def find_l1_bisector(p1, p2, width, height):
    """Generate L1 bisector between two sites"""
    x_distance = p1["site"][0] - p2["site"][0]
    y_distance = p1["site"][1] - p2["site"][1]

    mid_x = (p1["site"][0] + p2["site"][0]) / 2
    mid_y = (p1["site"][1] + p2["site"][1]) / 2

    vertexes = []
    up = None

    if same_point(p1["site"], p2["site"]):
        raise ValueError(
            f"Duplicate point: Points {p1['site']} and {p2['site']} are duplicates. Please remove one."
        )

    if abs(x_distance) == 0:
        bisector = create_bisector([p1, p2], False)
        bisector["points"] = [[0, mid_y], [width, mid_y]]
        return bisector

    if abs(y_distance) == 0:
        bisector = create_bisector([p1, p2], True)
        bisector["points"] = [[mid_x, 0], [mid_x, height]]
        return bisector

    slope = -1 if y_distance / x_distance > 0 else 1
    intercept = mid_y - mid_x * slope

    if abs(x_distance) > abs(y_distance):
        v1 = [(p1["site"][1] - intercept) / slope, p1["site"][1]]
        v2 = [(p2["site"][1] - intercept) / slope, p2["site"][1]]
        vertexes = [v1, v2]
        up = True
    elif abs(x_distance) < abs(y_distance):
        v1 = [p1["site"][0], p1["site"][0] * slope + intercept]
        v2 = [p2["site"][0], p2["site"][0] * slope + intercept]
        vertexes = [v1, v2]
        up = False
    else:
        if slope == 1:
            v1 = [p1["site"][1] - intercept, p1["site"][1]]
            v2 = [p2["site"][1] - intercept, p2["site"][1]]
            vertexes = [v1, v2]
            up = True
        else:
            v1 = [p1["site"][0], -p1["site"][0] + intercept]
            v2 = [p2["site"][0], -p2["site"][0] + intercept]
            vertexes = [v1, v2]
            up = False

    bisector = create_bisector([p1, p2], up)

    if up:
        sorted_verts = sorted(vertexes, key=lambda p: p[1])
        bisector["points"] = sorted(
            [[sorted_verts[0][0], 0]] + sorted_verts + [[sorted_verts[1][0], height]],
            key=lambda p: p[1],
        )
    else:
        sorted_verts = sorted(vertexes, key=lambda p: p[0])
        bisector["points"] = sorted(
            [[0, sorted_verts[0][1]]] + sorted_verts + [[width, sorted_verts[1][1]]],
            key=lambda p: p[0],
        )

    return bisector


def is_new_bisector_upward(hop_to, hop_from, site, go_up):
    """Check if bisector is traveling upward or downward"""
    # if hop_to["site"][0] == hop_from["site"][0]:
    #     return site["site"][1] > hop_to["site"][1]

    if hop_to["site"][0] - site["site"][0] == 0:
        return site["site"][1] > hop_to["site"][1]

    slope = (hop_to["site"][1] - site["site"][1]) / (
        hop_to["site"][0] - site["site"][0]
    )
    intercept = hop_to["site"][1] - (slope * hop_to["site"][0])

    is_above_line = hop_from["site"][1] > (slope * hop_from["site"][0]) + intercept

    return is_above_line
