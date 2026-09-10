"""
Pure point/segment geometry shared across the Voronoi pipeline.
"""

import math


def angle(p1, p2):
    """Calculate angle between two points"""
    ang = math.atan2(p2[1] - p1[1], p2[0] - p1[0])
    if ang < 0:
        ang = math.pi + math.pi + ang
    return ang


def distance(p1, p2):
    """L1 distance between two points"""
    return abs(p1[0] - p2[0]) + abs(p1[1] - p2[1])


def same_point(p1, p2):
    """Check if two points are the same"""
    return p1[0] == p2[0] and p1[1] == p2[1]


def segment_intersection(l1, l2):
    """Find intersection of two line segments"""
    x0, y0 = l1[0]
    x1, y1 = l1[1]
    x2, y2 = l2[0]
    x3, y3 = l2[1]

    denom = (y3 - y2) * (x1 - x0) - (x3 - x2) * (y1 - y0)

    if denom == 0:
        return None

    ua = ((x3 - x2) * (y0 - y2) - (y3 - y2) * (x0 - x2)) / denom
    ub = ((x1 - x0) * (y0 - y2) - (y1 - y0) * (x0 - x2)) / denom

    if not (ua >= 0 and ua <= 1 and ub >= 0 and ub <= 1):
        return False

    return [x0 + ua * (x1 - x0), y0 + ua * (y1 - y0)]
