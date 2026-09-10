"""
Library facade. Re-exports the public API from the pipeline modules so
consumers keep importing from a single entry point.
"""

from geometry import distance, same_point, angle, segment_intersection
from l1_metric import find_l1_bisector
from preprocess import clean_data
from generator import generate_l1_voronoi

__all__ = [
    "generate_l1_voronoi",
    "clean_data",
    "distance",
    "same_point",
    "angle",
    "segment_intersection",
    "find_l1_bisector",
]
