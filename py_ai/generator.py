"""
Public generator entry points. generate_l1_voronoi composes the pipeline
stages (preprocess -> sort/init -> divide & conquer -> polygonize); the
exported functions are the library facade.
"""

from divide_conquer import recursive_split
from polygonizer import polygonize_site
from cell import to_cell
from preprocess import clean_data
from l1_metric import create_l1_metric
from bisector import create_site


def generate_l1_voronoi(site_points, width, height, nudge_data=True):
    """Generate an L1 Voronoi diagram"""
    if nudge_data:
        site_points = clean_data(site_points[:])

    sorted_points = sorted(site_points, key=lambda p: (p[0], p[1]))
    sites = [create_site(e) for e in sorted_points]

    metric = create_l1_metric(width, height)
    graph = recursive_split(sites, metric)

    polygonized = [polygonize_site(site, metric) for site in graph]

    return [to_cell(site) for site in polygonized]
