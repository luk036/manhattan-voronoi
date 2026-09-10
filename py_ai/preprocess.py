"""
Input preprocessing: nudges degenerate point configurations (duplicate
points, points on a square) that the algorithm cannot handle exactly.
"""


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
