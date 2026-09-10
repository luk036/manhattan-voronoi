"""
Presentation boundary: map the internal site/bisector graph onto plain result
objects. Consumers and adapters (SVG, JSON, Canvas) see these DTOs, so the
live, mutable graph stays an implementation detail.
"""


def to_svg_path(polygon_points):
    """Render an SVG path string for a polygon point list."""
    return "M " + " L".join(" ".join(str(c) for c in p) for p in polygon_points) + " Z"


def to_bisector_dto(bisector):
    """Map an internal bisector onto a plain result object."""
    return {
        "sites": [{"site": site["site"][:]} for site in bisector["sites"]],
        "up": bisector["up"],
        "points": [point[:] for point in bisector["points"]],
        "intersections": [point[:] for point in bisector["intersections"]],
        "compound": bisector["compound"],
        "merge_line": bisector["merge_line"] if "merge_line" in bisector else None,
    }


def to_cell(site):
    """Map an internal site onto its public cell, applying the SVG adapter."""
    cell = {
        "site": site["site"],
        "bisectors": [to_bisector_dto(b) for b in site["bisectors"]],
    }

    if "polygon_points" in site:
        cell["polygon_points"] = [point[:] for point in site["polygon_points"]]
        cell["d"] = to_svg_path(site["polygon_points"])
        cell["neighbors"] = [neighbor[:] for neighbor in site["neighbors"]]

    return cell
