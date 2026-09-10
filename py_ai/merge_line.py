"""
The merge-line traversal of the Lee & Wong divide-and-conquer merge. Walks a
bisector between two merged halves, hopping across the nearest intersecting
bisector on either side until the merge line exits the canvas (or an orphaned
bisector needs to be unwound).
"""

from geometry import angle, same_point
from bisector import (
    bisector_intersection,
    trim_bisector,
    find_hop_to,
    is_bisector_trapped,
    get_extreme_point,
    remove_bisector,
)


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


def walk_merge_line(
    current_r,
    current_l,
    current_bisector,
    current_crop_point,
    go_up,
    metric,
    crossed_bisector=None,
    merge_array=None,
):
    """
    Walk along merge line to find intersections

    :param current_r: Site
    :param current_l: Site
    :param current_bisector: Bisector
    :param current_crop_point: [x, y]
    :param go_up: boolean
    :param metric: Metric strategy object
    :param crossed_bisector: Bisector or None
    :param merge_array: list of Bisectors
    :returns: list of Bisectors
    """
    if merge_array is None:
        merge_array = []

    # Ensure both sites are in current bisector
    if not all(e == current_r or e == current_l for e in current_bisector["sites"]):
        current_bisector = metric["bisector"](current_r, current_l)
        trim_bisector(current_bisector, crossed_bisector, current_crop_point, metric)
        merge_array.append(current_bisector)

    # Process left bisectors
    crop_l_array = []
    for e in current_l["bisectors"]:
        point = bisector_intersection(current_bisector, e)
        if point:
            hop_to = next((d for d in e["sites"] if d != current_l), None)
            if hop_to and go_up == metric["is_upward"](
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
        new_merge_line = metric["bisector"](current_r, hop_to)
        trim_bisector(new_merge_line, e["bisector"], e["point"], metric)

        # Check if trapped
        is_trapped = all(
            not is_bisector_trapped(
                find_hop_to(d["bisector"], current_l), new_merge_line, metric
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
            if hop_to and go_up == metric["is_upward"](
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
        new_merge_line = metric["bisector"](current_l, hop_to)
        trim_bisector(new_merge_line, e["bisector"], e["point"], metric)

        is_trapped = all(
            not is_bisector_trapped(
                find_hop_to(d["bisector"], current_r), new_merge_line, metric
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
        left_orphan = check_for_orphans(current_r, current_l, go_up, metric)
        right_orphan = check_for_orphans(current_l, current_r, go_up, metric)

        if left_orphan:
            remove_bisector(left_orphan)

            hop_to = find_hop_to(left_orphan, current_l)
            current_r = find_correct_w(current_r, hop_to, metric)
            new_merge_bisector = metric["bisector"](hop_to, current_r)
            merge_array.append(new_merge_bisector)

            return walk_merge_line(
                current_r,
                hop_to,
                new_merge_bisector,
                current_crop_point,
                go_up,
                metric,
                crossed_bisector,
                merge_array,
            )
        elif right_orphan:
            remove_bisector(right_orphan)

            hop_to = find_hop_to(right_orphan, current_r)
            current_l = find_correct_w(current_l, hop_to, metric)
            new_merge_bisector = metric["bisector"](hop_to, current_l)
            merge_array.append(new_merge_bisector)

            return walk_merge_line(
                hop_to,
                current_l,
                new_merge_bisector,
                current_crop_point,
                go_up,
                metric,
                crossed_bisector,
                merge_array,
            )

        return merge_array

    # Determine which point to cross first
    cross = determine_first_border_cross(crop_r, crop_l, current_crop_point)

    if cross == "right":
        trim_bisector(crop_r["bisector"], current_bisector, crop_r["point"], metric)
        trim_bisector(current_bisector, crop_r["bisector"], crop_r["point"], metric)
        current_bisector["intersections"].append(crop_r["point"])
        crossed_bisector = crop_r["bisector"]
        current_r = next(
            (s for s in crop_r["bisector"]["sites"] if s != current_r), current_r
        )
        current_crop_point = crop_r["point"]
    elif cross == "left":
        trim_bisector(crop_l["bisector"], current_bisector, crop_l["point"], metric)
        trim_bisector(current_bisector, crop_l["bisector"], crop_l["point"], metric)
        current_bisector["intersections"].append(crop_l["point"])
        crossed_bisector = crop_l["bisector"]
        current_l = next(
            (s for s in crop_l["bisector"]["sites"] if s != current_l), current_l
        )
        current_crop_point = crop_l["point"]
    else:
        # Both
        if crop_r["bisector"]:
            trim_bisector(crop_r["bisector"], current_bisector, crop_r["point"], metric)
            trim_bisector(current_bisector, crop_r["bisector"], crop_r["point"], metric)
            current_bisector["intersections"].append(crop_r["point"])
            crossed_bisector = crop_r["bisector"]
            current_r = next(
                (s for s in crop_r["bisector"]["sites"] if s != current_r), current_r
            )
            current_crop_point = crop_r["point"]

        if crop_l["bisector"]:
            trim_bisector(crop_l["bisector"], current_bisector, crop_l["point"], metric)
            trim_bisector(current_bisector, crop_l["bisector"], crop_l["point"], metric)
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
        metric,
        crossed_bisector,
        merge_array,
    )


def find_correct_w(w, nearest_neighbor, metric):
    """Ensure starting point is correct and would not result in trapped bisector"""
    starting_bisector = metric["bisector"](w, nearest_neighbor)

    w_traps = []
    for e in w["bisectors"]:
        hop_to = find_hop_to(e, w)
        is_trapped = is_bisector_trapped(hop_to, starting_bisector, metric)
        if is_trapped:
            w_traps.append({"hop_to": hop_to, "is_trapped": is_trapped})

    w_traps.sort(
        key=lambda x: metric["distance"](
            x["hop_to"]["site"], nearest_neighbor["site"]
        )
    )

    if w_traps:
        return find_correct_w(w_traps[0]["hop_to"], nearest_neighbor, metric)
    else:
        return w


def check_for_orphans(trapper, trapped, go_up, metric):
    """Check for orphaned bisectors"""
    orphans = []
    for bisector in trapped["bisectors"]:
        hop_to = find_hop_to(bisector, trapped)
        if go_up == (hop_to["site"][1] < trapped["site"][1]):
            if is_bisector_trapped(trapper, bisector, metric):
                orphans.append(bisector)

    if not orphans:
        return None

    orphans.sort(
        key=lambda a: _orphan_sort_key(a, trapped, trapper, go_up, metric)
    )
    return orphans[0] if orphans else None


def _orphan_sort_key(bisector, trapped, trapper, go_up, metric):
    """Helper for sorting orphans"""
    hop_to_a = find_hop_to(bisector, trapped)
    merge_line_a = metric["bisector"](hop_to_a, trapper)
    extreme_a = get_extreme_point(merge_line_a, go_up)

    def get_sort_val(b):
        ht = find_hop_to(b, trapped)
        ml = metric["bisector"](ht, trapper)
        return get_extreme_point(ml, go_up)

    # Simplified - just return extreme point
    return extreme_a


def determine_starting_bisector(w, nearest_neighbor, metric, last_intersect=None):
    """Determine starting bisector for the merge process"""
    z = [metric["width"], w["site"][1]]

    if not last_intersect:
        last_intersect = w["site"]

    zline = {"points": [w["site"], z]}

    intersection = None
    for bisector in nearest_neighbor["bisectors"]:
        pt = bisector_intersection(zline, bisector)
        if pt:
            intersection = {"point": pt, "bisector": bisector}
            break

    if intersection and metric["distance"](
        w["site"], intersection["point"]
    ) > metric["distance"](nearest_neighbor["site"], intersection["point"]):
        starting_bisector = metric["bisector"](w, nearest_neighbor)
        return {
            "starting_bisector": starting_bisector,
            "w": w,
            "nearest_neighbor": nearest_neighbor,
            "starting_intersection": intersection["point"]
            if intersection
            else w["site"],
        }
    elif intersection and metric["distance"](
        w["site"], intersection["point"]
    ) < metric["distance"](nearest_neighbor["site"], intersection["point"]):
        if intersection["point"][0] > last_intersect[0]:
            next_r = next(
                (e for e in intersection["bisector"]["sites"] if e != nearest_neighbor),
                nearest_neighbor,
            )
            return determine_starting_bisector(
                w, next_r, metric, intersection["point"]
            )

    w = find_correct_w(w, nearest_neighbor, metric)
    starting_bisector = metric["bisector"](w, nearest_neighbor)

    return {
        "starting_bisector": starting_bisector,
        "w": w,
        "nearest_neighbor": nearest_neighbor,
        "starting_intersection": intersection["point"] if intersection else w["site"],
    }
