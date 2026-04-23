"""
Generate L1 (Manhattan distance) Voronoi diagrams using Lee and Wong's algorithm.
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


def clean_data(data):
    """
    Nudge points to hopefully eliminate square bisectors

    :param data: list of points [[x,y], ...]
    :returns: modified data list
    """
    for i, e in enumerate(data):
        for j, d in enumerate(data):
            if i != j and abs(d[0] - e[0]) == abs(d[1] - e[1]):
                d[0] = d[0] + 1e-10 * d[1]
                d[1] = d[1] + 2e-10 * d[0]
    return data


def generate_l1_voronoi(site_points, width, height, nudge_data=True):
    """Generate an L1 Voronoi diagram"""
    if nudge_data:
        site_points = clean_data(site_points[:])

    sorted_points = sorted(site_points, key=lambda p: (p[0], p[1]))
    sites = [{"site": e, "bisectors": []} for e in sorted_points]

    find_bisector = curry_find_bisector(find_l1_bisector, width, height)
    graph = recursive_split(sites, find_bisector, width, height)

    graph = finalize_sites(graph, width, height)

    return graph


def recursive_split(split_array, find_bisector, width, height):
    """
    Recursively split and merge sets of points

    :param split_array: list of Site objects
    :param find_bisector: function
    :param width: integer
    :param height: integer
    :returns: list of Site objects
    """
    # If more than two points, split recursively
    if len(split_array) > 2:
        split_point = (len(split_array) - len(split_array) % 2) // 2

        # Merge the child diagrams
        left = recursive_split(split_array[:split_point], find_bisector, width, height)
        right = recursive_split(split_array[split_point:], find_bisector, width, height)

        # The current working sites
        right_sorted = sorted(
            right, key=lambda s: distance(left[-1]["site"], s["site"])
        )

        starting_info = determine_starting_bisector(
            left[-1], right_sorted[0], width, None, find_bisector
        )

        initial_bisector = starting_info["starting_bisector"]
        initial_r = starting_info["nearest_neighbor"]
        initial_l = starting_info["w"]

        up_stroke_array = walk_merge_line(
            initial_r,
            initial_l,
            initial_bisector,
            [width, height],
            True,
            None,
            [],
            find_bisector,
        )
        down_stroke_array = walk_merge_line(
            initial_r,
            initial_l,
            initial_bisector,
            [0, 0],
            False,
            None,
            [],
            find_bisector,
        )

        # Combine all merge arrays
        merge_array = [initial_bisector] + up_stroke_array + down_stroke_array

        for bisector in merge_array:
            bisector["merge_line"] = len(split_array)
            bisector["sites"][0]["bisectors"] = clear_out_orphans(
                bisector["sites"][0], bisector["sites"][1]
            )
            bisector["sites"][1]["bisectors"] = clear_out_orphans(
                bisector["sites"][1], bisector["sites"][0]
            )

            for site in bisector["sites"]:
                site["bisectors"].append(bisector)

        return left + right

    # Otherwise, determine the vertices if it has two sites
    elif len(split_array) == 2:
        bisector = find_bisector(split_array[0], split_array[1])
        for e in split_array:
            e["bisectors"].append(bisector)
        return split_array

    # If it has just one, just return it
    else:
        return split_array


def walk_merge_line(
    current_r,
    current_l,
    current_bisector,
    current_crop_point,
    go_up,
    crossed_bisector=None,
    merge_array=None,
    find_bisector=None,
):
    """
    Walk along merge line to find intersections

    :param current_r: Site
    :param current_l: Site
    :param current_bisector: Bisector
    :param current_crop_point: [x, y]
    :param go_up: boolean
    :param crossed_bisector: Bisector or None
    :param merge_array: list of Bisectors
    :param find_bisector: function
    :returns: list of Bisectors
    """
    if merge_array is None:
        merge_array = []

    # Ensure both sites are in current bisector
    if not all(e == current_r or e == current_l for e in current_bisector["sites"]):
        current_bisector = find_bisector(current_r, current_l)
        trim_bisector(current_bisector, crossed_bisector, current_crop_point)
        merge_array.append(current_bisector)

    # Process left bisectors
    crop_l_array = []
    for e in current_l["bisectors"]:
        point = bisector_intersection(current_bisector, e)
        if point:
            hop_to = next((d for d in e["sites"] if d != current_l), None)
            if hop_to and go_up == is_new_bisector_upward(
                hop_to, current_l, current_r, go_up
            ):
                if not same_point(point, current_crop_point) or e != crossed_bisector:
                    crop_l_array.append({"bisector": e, "point": point})

    # Sort by angle
    crop_l_array.sort(
        key=lambda x: angle(
            current_l["site"], find_hop_to(x["bisector"], current_l)["site"]
        )
    )

    # Filter trapped bisectors
    filtered_l = []
    for i, e in enumerate(crop_l_array):
        hop_to = find_hop_to(e["bisector"], current_l)
        new_merge_line = find_bisector(current_r, hop_to)
        trim_bisector(new_merge_line, e["bisector"], e["point"])

        # Check if trapped
        is_trapped = all(
            not is_bisector_trapped(
                find_hop_to(d["bisector"], current_l), new_merge_line
            )
            or find_hop_to(d["bisector"], current_l) == hop_to
            for d in crop_l_array[:i]
        )
        if is_trapped:
            filtered_l.append(e)
    crop_l_array = filtered_l

    # Process right bisectors
    crop_r_array = []
    for e in current_r["bisectors"]:
        point = bisector_intersection(current_bisector, e)
        if point:
            hop_to = next((d for d in e["sites"] if d != current_r), None)
            if hop_to and go_up == is_new_bisector_upward(
                hop_to, current_r, current_l, go_up
            ):
                if not same_point(point, current_crop_point) or e != crossed_bisector:
                    crop_r_array.append({"bisector": e, "point": point})

    crop_r_array.sort(
        key=lambda x: angle(
            current_r["site"], find_hop_to(x["bisector"], current_r)["site"]
        )
    )

    filtered_r = []
    for i, e in enumerate(crop_r_array):
        hop_to = find_hop_to(e["bisector"], current_r)
        new_merge_line = find_bisector(current_l, hop_to)
        trim_bisector(new_merge_line, e["bisector"], e["point"])

        is_trapped = all(
            not is_bisector_trapped(
                find_hop_to(d["bisector"], current_r), new_merge_line
            )
            or find_hop_to(d["bisector"], current_r) == hop_to
            for d in crop_r_array[:i]
        )
        if is_trapped:
            filtered_r.append(e)
    crop_r_array = filtered_r

    # Determine crop points
    infinity_pt = (
        [float("inf"), float("inf")] if go_up else [-float("inf"), -float("inf")]
    )
    crop_l = (
        crop_l_array[0]
        if crop_l_array and crop_l_array[0] != current_bisector
        else {"bisector": None, "point": infinity_pt}
    )
    crop_r = (
        crop_r_array[0]
        if crop_r_array and crop_r_array[0] != current_bisector
        else {"bisector": None, "point": infinity_pt}
    )

    # If no intersection, we're done
    if not crop_l["bisector"] and not crop_r["bisector"]:
        # Check for orphans
        left_orphan = check_for_orphans(current_r, current_l, go_up, find_bisector)
        right_orphan = check_for_orphans(current_l, current_r, go_up, find_bisector)

        if left_orphan:
            for site in left_orphan["sites"]:
                site["bisectors"] = [b for b in site["bisectors"] if b != left_orphan]

            hop_to = find_hop_to(left_orphan, current_l)
            current_r = find_correct_w(current_r, hop_to, find_bisector)
            new_merge_bisector = find_bisector(hop_to, current_r)
            merge_array.append(new_merge_bisector)

            return walk_merge_line(
                current_r,
                hop_to,
                new_merge_bisector,
                current_crop_point,
                go_up,
                crossed_bisector,
                merge_array,
                find_bisector,
            )
        elif right_orphan:
            for site in right_orphan["sites"]:
                site["bisectors"] = [b for b in site["bisectors"] if b != right_orphan]

            hop_to = find_hop_to(right_orphan, current_r)
            current_l = find_correct_w(current_l, hop_to, find_bisector)
            new_merge_bisector = find_bisector(hop_to, current_l)
            merge_array.append(new_merge_bisector)

            return walk_merge_line(
                hop_to,
                current_l,
                new_merge_bisector,
                current_crop_point,
                go_up,
                crossed_bisector,
                merge_array,
                find_bisector,
            )

        return merge_array

    # Determine which point to cross first
    cross = determine_first_border_cross(crop_r, crop_l, current_crop_point)

    if cross == "right":
        trim_bisector(crop_r["bisector"], current_bisector, crop_r["point"])
        trim_bisector(current_bisector, crop_r["bisector"], crop_r["point"])
        current_bisector["intersections"].append(crop_r["point"])
        crossed_bisector = crop_r["bisector"]
        current_r = next(
            (s for s in crop_r["bisector"]["sites"] if s != current_r), current_r
        )
        current_crop_point = crop_r["point"]
    elif cross == "left":
        trim_bisector(crop_l["bisector"], current_bisector, crop_l["point"])
        trim_bisector(current_bisector, crop_l["bisector"], crop_l["point"])
        current_bisector["intersections"].append(crop_l["point"])
        crossed_bisector = crop_l["bisector"]
        current_l = next(
            (s for s in crop_l["bisector"]["sites"] if s != current_l), current_l
        )
        current_crop_point = crop_l["point"]
    else:
        # Both
        if crop_r["bisector"]:
            trim_bisector(crop_r["bisector"], current_bisector, crop_r["point"])
            trim_bisector(current_bisector, crop_r["bisector"], crop_r["point"])
            current_bisector["intersections"].append(crop_r["point"])
            crossed_bisector = crop_r["bisector"]
            current_r = next(
                (s for s in crop_r["bisector"]["sites"] if s != current_r), current_r
            )
            current_crop_point = crop_r["point"]

        if crop_l["bisector"]:
            trim_bisector(crop_l["bisector"], current_bisector, crop_l["point"])
            trim_bisector(current_bisector, crop_l["bisector"], crop_l["point"])
            current_bisector["intersections"].append(crop_l["point"])
            crossed_bisector = crop_l["bisector"]
            current_l = next(
                (s for s in crop_l["bisector"]["sites"] if s != current_l), current_l
            )
            current_crop_point = crop_l["point"]

    return walk_merge_line(
        current_r,
        current_l,
        current_bisector,
        current_crop_point,
        go_up,
        crossed_bisector,
        merge_array,
        find_bisector,
    )


def find_correct_w(w, nearest_neighbor, find_bisector):
    """Ensure starting point is correct and would not result in trapped bisector"""
    starting_bisector = find_bisector(w, nearest_neighbor)

    w_traps = []
    for e in w["bisectors"]:
        hop_to = find_hop_to(e, w)
        is_trapped = is_bisector_trapped(hop_to, starting_bisector)
        if is_trapped:
            w_traps.append({"hop_to": hop_to, "is_trapped": is_trapped})

    w_traps.sort(key=lambda x: distance(x["hop_to"]["site"], nearest_neighbor["site"]))

    if w_traps:
        return find_correct_w(w_traps[0]["hop_to"], nearest_neighbor, find_bisector)
    else:
        return w


def check_for_orphans(trapper, trapped, go_up, find_bisector):
    """Check for orphaned bisectors"""
    orphans = []
    for bisector in trapped["bisectors"]:
        hop_to = find_hop_to(bisector, trapped)
        if go_up == (hop_to["site"][1] < trapped["site"][1]):
            if is_bisector_trapped(trapper, bisector):
                orphans.append(bisector)

    if not orphans:
        return None

    orphans.sort(
        key=lambda a: _orphan_sort_key(a, trapped, trapper, go_up, find_bisector)
    )
    return orphans[0] if orphans else None


def _orphan_sort_key(bisector, trapped, trapper, go_up, find_bisector):
    """Helper for sorting orphans"""
    hop_to_a = find_hop_to(bisector, trapped)
    merge_line_a = find_bisector(hop_to_a, trapper)
    extreme_a = get_extreme_point(merge_line_a, go_up)

    def get_sort_val(b):
        ht = find_hop_to(b, trapped)
        ml = find_bisector(ht, trapper)
        return get_extreme_point(ml, go_up)

    # Simplified - just return extreme point
    return extreme_a


def determine_first_border_cross(crop_r, crop_l, current_crop_point):
    """Determine which border to cross first"""
    if abs(crop_r["point"][1] - current_crop_point[1]) == abs(
        crop_l["point"][1] - current_crop_point[1]
    ):
        return None
    elif abs(crop_r["point"][1] - current_crop_point[1]) < abs(
        crop_l["point"][1] - current_crop_point[1]
    ):
        return "right"
    else:
        return "left"


def determine_starting_bisector(
    w, nearest_neighbor, width, last_intersect=None, find_bisector=None
):
    """Determine starting bisector for the merge process"""
    z = [width, w["site"][1]]

    if not last_intersect:
        last_intersect = w["site"]

    zline = {"points": [w["site"], z]}

    intersection = None
    for bisector in nearest_neighbor["bisectors"]:
        pt = bisector_intersection(zline, bisector)
        if pt:
            intersection = {"point": pt, "bisector": bisector}
            break

    if intersection and distance(w["site"], intersection["point"]) > distance(
        nearest_neighbor["site"], intersection["point"]
    ):
        starting_bisector = find_bisector(w, nearest_neighbor)
        return {
            "starting_bisector": starting_bisector,
            "w": w,
            "nearest_neighbor": nearest_neighbor,
            "starting_intersection": intersection["point"]
            if intersection
            else w["site"],
        }
    elif intersection and distance(w["site"], intersection["point"]) < distance(
        nearest_neighbor["site"], intersection["point"]
    ):
        if intersection["point"][0] > last_intersect[0]:
            next_r = next(
                (e for e in intersection["bisector"]["sites"] if e != nearest_neighbor),
                nearest_neighbor,
            )
            return determine_starting_bisector(
                w, next_r, width, intersection["point"], find_bisector
            )

    w = find_correct_w(w, nearest_neighbor, find_bisector)
    starting_bisector = find_bisector(w, nearest_neighbor)

    return {
        "starting_bisector": starting_bisector,
        "w": w,
        "nearest_neighbor": nearest_neighbor,
        "starting_intersection": intersection["point"] if intersection else w["site"],
    }


# === Helper Functions ===


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


def find_hop_to(bisector, hop_from):
    """Find the other point across a bisector"""
    for e in bisector["sites"]:
        if e != hop_from:
            return e
    return bisector["sites"][0]


def curry_find_bisector(callback, width, height):
    """Curry find bisector function with width and height"""

    def curried(p1, p2):
        return callback(p1, p2, width, height)

    return curried


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
        vertexes = [[0, mid_y], [width, mid_y]]
        return {
            "sites": [p1, p2],
            "up": False,
            "points": vertexes,
            "intersections": [],
            "compound": False,
        }

    if abs(y_distance) == 0:
        vertexes = [[mid_x, 0], [mid_x, height]]
        return {
            "sites": [p1, p2],
            "up": True,
            "points": vertexes,
            "intersections": [],
            "compound": False,
        }

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

    bisector = {
        "sites": [p1, p2],
        "up": up,
        "points": [],
        "intersections": [],
        "compound": False,
    }

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


def clear_out_orphans(orphanage, trap_point):
    """Clear out orphans when a new merge line is created"""
    return [b for b in orphanage["bisectors"] if not is_bisector_trapped(trap_point, b)]


def is_bisector_trapped(trap_point, bisector):
    """Determine if bisector is trapped in a site's polygon"""
    site0 = bisector["sites"][0]
    site1 = bisector["sites"][1]

    return all(
        distance(trap_point["site"], point) <= distance(site0["site"], point)
        and distance(trap_point["site"], point) <= distance(site1["site"], point)
        for point in bisector["points"]
    )


def get_extreme_point(bisector, go_up):
    """Find highest or lowest point of bisector"""
    if go_up:
        return max(point[1] for point in bisector["points"])
    else:
        return min(point[1] for point in bisector["points"])


def trim_bisector(target, intersector, intersection):
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
        if distance(point, target["sites"][0]["site"])
        < distance(point, polygon_site["site"])
        and distance(point, target["sites"][1]["site"])
        < distance(point, polygon_site["site"])
    ]
    new_points.append(intersection)

    # Sort
    if target["up"]:
        target["points"] = sorted(new_points, key=lambda p: p[1])
    else:
        target["points"] = sorted(new_points, key=lambda p: p[0])


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


# === Post-processing to add polygon points and SVG d ===


def finalize_sites(sites, width, height):
    """Add polygonPoints, d, and neighbors to sites"""
    corners = [[0, 0], [width, 0], [width, height], [0, height]]

    result = []
    for site in sites:
        # Build polygon from bisectors
        bisectors = site["bisectors"]

        if not bisectors:
            result.append(site)
            continue

        # Find starting bisector
        start_bisector = None
        for b in bisectors:
            if any(is_point_on_edge(p, width, height) for p in b["points"]):
                start_bisector = b
                break
        start_bisector = start_bisector or bisectors[0]

        starting_points = start_bisector["points"]

        # Reverse if ends on edge
        if is_point_on_edge(starting_points[-1], width, height):
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

                d1 = distance(last, b["points"][0])
                d2 = distance(last, b["points"][-1])
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

        # Handle case where polygon starts and ends on different edges
        if is_point_on_edge(polygon_points[0], width, height) and is_point_on_edge(
            polygon_points[-1], width, height
        ):
            if not _are_points_on_same_edge(
                polygon_points[0], polygon_points[-1], width, height
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

        # Sort by angle
        polygon_points = sorted(polygon_points, key=lambda p: angle(site["site"], p))

        # Create SVG d
        d_str = (
            "M " + " L".join(" ".join(str(c) for c in p) for p in polygon_points) + " Z"
        )

        # Get neighbors
        neighbors = [find_hop_to(b, site["site"]) for b in site["bisectors"]]

        site["polygon_points"] = polygon_points
        site["d"] = d_str
        site["neighbors"] = neighbors
        result.append(site)

    return result


def is_point_on_edge(point, width, height):
    """Check if point is on an edge"""
    return point[0] == 0 or point[0] == width or point[1] == 0 or point[1] == height


def _are_points_on_same_edge(p1, p2, width, height):
    """Check if two points are on the same edge"""
    return (
        (p1[0] == p2[0] and p1[0] == 0)
        or (p1[0] == p2[0] and p1[0] == width)
        or (p1[1] == p2[1] and p1[1] == 0)
        or (p1[1] == p2[1] and p1[1] == height)
    )


# === Exports ===

__all__ = [
    "generate_voronoi_points",
    "generate_l1_voronoi",
    "clean_data",
    "finalize_sites",
]
