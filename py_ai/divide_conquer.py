"""
The divide step of the algorithm: recursively split a sorted site list into
singletons, building base bisectors for pairs and walking merge lines to
combine sibling halves.
"""

from bisector import link_bisector_to_sites, clear_out_orphans
from merge_line import walk_merge_line, determine_starting_bisector


def recursive_split(split_array, metric):
    """
    Recursively split and merge sets of points

    :param split_array: list of Site objects
    :param metric: Metric strategy object
    :returns: list of Site objects
    """
    # If more than two points, split recursively
    if len(split_array) > 2:
        split_point = (len(split_array) - len(split_array) % 2) // 2

        # Merge the child diagrams
        left = recursive_split(split_array[:split_point], metric)
        right = recursive_split(split_array[split_point:], metric)

        # The current working sites
        right_sorted = sorted(
            right, key=lambda s: metric["distance"](left[-1]["site"], s["site"])
        )

        starting_info = determine_starting_bisector(left[-1], right_sorted[0], metric)

        initial_bisector = starting_info["starting_bisector"]
        initial_r = starting_info["nearest_neighbor"]
        initial_l = starting_info["w"]

        up_stroke_array = walk_merge_line(
            initial_r,
            initial_l,
            initial_bisector,
            [metric["width"], metric["height"]],
            True,
            metric,
        )
        down_stroke_array = walk_merge_line(
            initial_r,
            initial_l,
            initial_bisector,
            [0, 0],
            False,
            metric,
        )

        # Combine all merge arrays
        merge_array = [initial_bisector] + up_stroke_array + down_stroke_array

        for bisector in merge_array:
            bisector["merge_line"] = len(split_array)
            bisector["sites"][0]["bisectors"] = clear_out_orphans(
                bisector["sites"][0], bisector["sites"][1], metric
            )
            bisector["sites"][1]["bisectors"] = clear_out_orphans(
                bisector["sites"][1], bisector["sites"][0], metric
            )
            link_bisector_to_sites(bisector)

        return left + right

    # Otherwise, determine the vertices if it has two sites
    elif len(split_array) == 2:
        bisector = metric["bisector"](split_array[0], split_array[1])
        link_bisector_to_sites(bisector)
        return split_array

    # If it has just one, just return it
    else:
        return split_array
