import math
import random
from typing import List, Tuple, Optional, Dict, Any
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.path import Path
import numpy as np


class Site:
    def __init__(self, site: List[float]):
        self.site = site
        self.bisectors = []
        self.polygon_points = []
        self.d = ""
        self.neighbors = []


class Bisector:
    def __init__(self, sites: List[Site], up: bool, points: List[List[float]], 
                 intersections: List[List[float]], compound: bool = False):
        self.sites = sites
        self.up = up
        self.points = points
        self.intersections = intersections
        self.compound = compound
        self.merge_line = None


def generate_voronoi_points(points: List[List[float]], width: int, height: int, 
                           distance_callback):
    """Generate Voronoi points via a basic, naive algorithm"""
    colors = [{'point': p, 'color': [random.randint(0, 255) for _ in range(3)]} 
              for p in points]
    
    image_data = []
    for x in range(width):
        for y in range(height):
            coordinate = [x, y]
            closest = colors[0]
            
            for color in colors[1:]:
                dist_closest = distance_callback(closest['point'], coordinate)
                dist_color = distance_callback(color['point'], coordinate)
                
                if dist_color < dist_closest:
                    closest = color
                elif dist_color == dist_closest:
                    # Tie case - set to black
                    closest = {'color': [0, 0, 0]}
            
            image_data.append(closest['color'])
    
    return image_data


def clean_data(data: List[List[float]]):
    """Nudge points to hopefully eliminate square bisectors"""
    for i, e in enumerate(data):
        for j, d in enumerate(data):
            if i != j and abs(d[0] - e[0]) == abs(d[1] - e[1]):
                d[0] = d[0] + 1e-10 * d[1]
                d[1] = d[1] + 2e-10 * d[0]
    return data


def generate_l1_voronoi(site_points: List[List[float]], width: int, height: int, 
                       nudge_data: bool = True):
    """Generate an L1 Voronoi diagram"""
    if nudge_data:
        site_points = clean_data(site_points)
    
    # Sort points by x axis, breaking ties with y
    sites = [{'site': p, 'bisectors': []} for p in site_points]
    sites.sort(key=lambda p: (p['site'][0], p['site'][1]))
    
    # Convert to Site objects
    site_objects = [Site(s['site']) for s in sites]
    
    find_bisector = curry_find_bisector(find_l1_bisector, width, height)
    graph = recursive_split(site_objects, find_bisector, width, height)
    
    for site in graph:
        # Build polygon points from bisectors
        if not site.bisectors:
            continue
            
        # Find a bisector on an edge if we have one
        start_bisector = next((b for b in site.bisectors 
                              if any(is_point_on_edge(p, width, height) for p in b.points)), 
                             site.bisectors[0])
        
        starting_points = start_bisector.points.copy()
        
        if is_point_on_edge(starting_points[-1], width, height):
            starting_points.reverse()
        
        total = {'points': starting_points, 'used': [start_bisector]}
        
        # Process remaining bisectors
        for bisector in site.bisectors[1:]:
            if bisector in total['used']:
                continue
                
            last = total['points'][-1]
            
            # Find closest connection point
            next_bisector = min((b for b in site.bisectors if b not in total['used']), 
                               key=lambda b: min(manhattan_distance(last, b.points[0]), 
                                               manhattan_distance(last, b.points[-1])))
            
            next_points = next_bisector.points.copy()
            
            if same_point(next_points[-1], last):
                next_points.reverse()
            
            total['points'].extend(next_points)
            total['used'].append(next_bisector)
        
        site.polygon_points = total['points']
        
        # Handle edge cases
        corners = [[0, 0], [width, 0], [width, height], [0, height]]
        
        if (is_point_on_edge(site.polygon_points[0], width, height) and 
            is_point_on_edge(site.polygon_points[-1], width, height) and
            not are_points_on_same_edge(site.polygon_points[0], 
                                       site.polygon_points[-1], width, height)):
            
            # Create temporary bisector for corner check
            temp_corners = []
            for c in corners:
                temp_site = Site(c)
                temp_bisector = find_bisector(site, temp_site)
                if all(not bisector_intersection(temp_bisector, b) for b in site.bisectors):
                    temp_corners.append(c)
            
            filtered_corners = temp_corners
            
            site.polygon_points.extend(filtered_corners)
        
        # Sort points by angle around the site
        site.polygon_points.sort(key=lambda p: angle(site.site, p))
        
        # Create SVG path string
        site.d = f"M {' '.join(f'{x} {y}' for x, y in site.polygon_points)} Z"
        
        # Set neighbors
        site.neighbors = [find_hop_to(b, site).site for b in site.bisectors]
    
    return graph


def is_point_on_edge(point: List[float], width: int, height: int) -> bool:
    return (point[0] == 0 or point[0] == width or 
            point[1] == 0 or point[1] == height)


def are_points_on_same_edge(p1: List[float], p2: List[float], width: int, height: int) -> bool:
    return ((p1[0] == p2[0] and p1[0] == 0) or
            (p1[0] == p2[0] and p1[0] == width) or
            (p1[1] == p2[1] and p1[1] == 0) or
            (p1[1] == p2[1] and p1[1] == height))


def recursive_split(split_array: List[Site], find_bisector, width: int, height: int):
    """Recursively split and merge sets of points"""
    if len(split_array) > 2:
        split_point = len(split_array) // 2
        
        L = recursive_split(split_array[:split_point], find_bisector, width, height)
        R = recursive_split(split_array[split_point:], find_bisector, width, height)
        
        # Find nearest neighbors for merging
        neighbor_array = sorted(R, key=lambda a: manhattan_distance(L[-1].site, a.site))
        
        starting_info = determine_starting_bisector(L[-1], neighbor_array[0], width, 
                                                   None, find_bisector)
        
        initial_bisector = starting_info['starting_bisector']
        initial_r = starting_info['nearest_neighbor']
        initial_l = starting_info['w']
        
        up_stroke_array = walk_merge_line(initial_r, initial_l, initial_bisector, 
                                         [width, height], True, None, [], find_bisector)
        down_stroke_array = walk_merge_line(initial_r, initial_l, initial_bisector, 
                                           [0, 0], False, None, [], find_bisector)
        
        merge_array = [initial_bisector] + up_stroke_array + down_stroke_array
        
        for bisector in merge_array:
            bisector.merge_line = len(split_array)
            
            bisector.sites[0].bisectors = clear_out_orphans(bisector.sites[0], 
                                                           bisector.sites[1])
            bisector.sites[1].bisectors = clear_out_orphans(bisector.sites[1], 
                                                           bisector.sites[0])
            
            for site in bisector.sites:
                site.bisectors.append(bisector)
        
        return L + R
    
    elif len(split_array) == 2:
        bisector = find_bisector(split_array[0], split_array[1])
        split_array[0].bisectors.append(bisector)
        split_array[1].bisectors.append(bisector)
        return split_array
    
    else:
        return split_array


def walk_merge_line(current_r: Site, current_l: Site, current_bisector: Bisector, 
                   current_crop_point: List[float], go_up: bool, crossed_border: Optional[Bisector], 
                   merge_array: List[Bisector], find_bisector):
    """Walk the merge line between two sets of sites"""
    if not all(site in current_bisector.sites for site in [current_r, current_l]):
        current_bisector = find_bisector(current_r, current_l)
        trim_bisector(current_bisector, crossed_border, current_crop_point)
        merge_array.append(current_bisector)
    
    # Process left side
    crop_l_array = []
    for bisector in current_l.bisectors:
        hop_to = find_hop_to(bisector, current_l)
        intersection = bisector_intersection(current_bisector, bisector)
        
        if (intersection and 
            go_up == is_new_bisector_upward(hop_to, current_l, current_r, go_up) and
            (not same_point(intersection, current_crop_point) or bisector != crossed_border)):
            crop_l_array.append({'bisector': bisector, 'point': intersection})
    
    crop_l_array.sort(key=lambda x: angle(current_l.site, find_hop_to(x['bisector'], current_l).site))
    
    # Filter trapped bisectors
    filtered_l_array = []
    for item in crop_l_array:
        hop_to = find_hop_to(item['bisector'], current_l)
        new_merge_line = find_bisector(current_r, hop_to)
        trim_bisector(new_merge_line, item['bisector'], item['point'])
        
        if all(not is_bisector_trapped(find_hop_to(other['bisector'], current_l), new_merge_line) 
               or find_hop_to(other['bisector'], current_l) == hop_to 
               for other in crop_l_array):
            filtered_l_array.append(item)
    
    # Process right side (similar to left)
    crop_r_array = []
    for bisector in current_r.bisectors:
        hop_to = find_hop_to(bisector, current_r)
        intersection = bisector_intersection(current_bisector, bisector)
        
        if (intersection and 
            go_up == is_new_bisector_upward(hop_to, current_r, current_l, go_up) and
            (not same_point(intersection, current_crop_point) or bisector != crossed_border)):
            crop_r_array.append({'bisector': bisector, 'point': intersection})
    
    crop_r_array.sort(key=lambda x: angle(current_r.site, find_hop_to(x['bisector'], current_r).site))
    
    filtered_r_array = []
    for item in crop_r_array:
        hop_to = find_hop_to(item['bisector'], current_r)
        new_merge_line = find_bisector(current_l, hop_to)
        trim_bisector(new_merge_line, item['bisector'], item['point'])
        
        if all(not is_bisector_trapped(find_hop_to(other['bisector'], current_r), new_merge_line) 
               or find_hop_to(other['bisector'], current_r) == hop_to 
               for other in crop_r_array):
            filtered_r_array.append(item)
    
    crop_l = filtered_l_array[0] if filtered_l_array and filtered_l_array[0]['bisector'] != current_bisector else {'bisector': None, 'point': [float('inf'), float('inf')] if go_up else [float('-inf'), float('-inf')]}
    crop_r = filtered_r_array[0] if filtered_r_array and filtered_r_array[0]['bisector'] != current_bisector else {'bisector': None, 'point': [float('inf'), float('inf')] if go_up else [float('-inf'), float('-inf')]}
    
    # If no intersection, we're done
    if not crop_l['bisector'] and not crop_r['bisector']:
        left_orphan = check_for_orphans(current_r, current_l, go_up, find_bisector)
        right_orphan = check_for_orphans(current_l, current_r, go_up, find_bisector)
        
        if left_orphan:
            # Remove trapped bisector
            for site in left_orphan.sites:
                site.bisectors = [b for b in site.bisectors if b != left_orphan]
            
            hop_to = find_hop_to(left_orphan, current_l)
            current_r = find_correct_w(current_r, hop_to, find_bisector)
            new_merge_bisector = find_bisector(hop_to, current_r)
            merge_array.append(new_merge_bisector)
            
            return walk_merge_line(current_r, hop_to, new_merge_bisector, current_crop_point, 
                                 go_up, crossed_border, merge_array, find_bisector)
        
        elif right_orphan:
            # Remove trapped bisector
            for site in right_orphan.sites:
                site.bisectors = [b for b in site.bisectors if b != right_orphan]
            
            hop_to = find_hop_to(right_orphan, current_r)
            current_l = find_correct_w(current_l, hop_to, find_bisector)
            new_merge_bisector = find_bisector(hop_to, current_l)
            merge_array.append(new_merge_bisector)
            
            return walk_merge_line(hop_to, current_l, new_merge_bisector, current_crop_point, 
                                 go_up, crossed_border, merge_array, find_bisector)
        
        return merge_array
    
    # Determine which point comes first
    first_cross = determine_first_border_cross(crop_r, crop_l, current_crop_point)
    
    if first_cross == "right":
        trim_bisector(crop_r['bisector'], current_bisector, crop_r['point'])
        trim_bisector(current_bisector, crop_r['bisector'], crop_r['point'])
        current_bisector.intersections.append(crop_r['point'])
        crossed_border = crop_r['bisector']
        current_r = next(s for s in crop_r['bisector'].sites if s != current_r)
        current_crop_point = crop_r['point']
    
    elif first_cross == "left":
        trim_bisector(crop_l['bisector'], current_bisector, crop_l['point'])
        trim_bisector(current_bisector, crop_l['bisector'], crop_l['point'])
        current_bisector.intersections.append(crop_l['point'])
        crossed_border = crop_l['bisector']
        current_l = next(s for s in crop_l['bisector'].sites if s != current_l)
        current_crop_point = crop_l['point']
    
    else:  # Both at same point
        trim_bisector(crop_r['bisector'], current_bisector, crop_r['point'])
        trim_bisector(current_bisector, crop_r['bisector'], crop_r['point'])
        current_bisector.intersections.append(crop_r['point'])
        crossed_border = crop_r['bisector']
        current_r = next(s for s in crop_r['bisector'].sites if s != current_r)
        current_crop_point = crop_r['point']
        
        trim_bisector(crop_l['bisector'], current_bisector, crop_l['point'])
        trim_bisector(current_bisector, crop_l['bisector'], crop_l['point'])
        current_bisector.intersections.append(crop_l['point'])
        crossed_border = crop_l['bisector']
        current_l = next(s for s in crop_l['bisector'].sites if s != current_l)
        current_crop_point = crop_l['point']
    
    return walk_merge_line(current_r, current_l, current_bisector, current_crop_point, 
                          go_up, crossed_border, merge_array, find_bisector)


def angle(p1: List[float], p2: List[float]) -> float:
    """Calculate angle between two points"""
    angle = math.atan2(p2[1] - p1[1], p2[0] - p1[0])
    
    if angle < 0:
        angle = math.pi + math.pi + angle
    
    return angle


def determine_first_border_cross(crop_r, crop_l, current_crop_point):
    """Determine which border crossing comes first"""
    r_dist = abs(crop_r['point'][1] - current_crop_point[1])
    l_dist = abs(crop_l['point'][1] - current_crop_point[1])
    
    if r_dist == l_dist:
        return None
    return "right" if r_dist < l_dist else "left"


def determine_starting_bisector(w: Site, nearest_neighbor: Site, width: int, 
                               last_intersect: Optional[List[float]], find_bisector):
    """Determine starting bisector for the merge process"""
    z = [width, w.site[1]]
    
    if not last_intersect:
        last_intersect = w.site
    
    # Create a temporary bisector for the zline
    temp_site = Site(z)
    zline = find_bisector(w, temp_site)
    
    intersection = None
    for bisector in nearest_neighbor.bisectors:
        point = bisector_intersection(zline, bisector)
        if point:
            intersection = {'point': point, 'bisector': bisector}
            break
    
    if (intersection and 
        manhattan_distance(w.site, intersection['point']) > 
        manhattan_distance(nearest_neighbor.site, intersection['point'])):
        
        starting_bisector = find_bisector(w, nearest_neighbor)
        return {
            'starting_bisector': starting_bisector,
            'w': w,
            'nearest_neighbor': nearest_neighbor,
            'starting_intersection': intersection['point'] if intersection['point'] else w.site
        }
    
    elif (intersection and 
          manhattan_distance(w.site, intersection['point']) < 
          manhattan_distance(nearest_neighbor.site, intersection['point']) and
          intersection['point'][0] > last_intersect[0]):
        
        next_r = next(s for s in intersection['bisector'].sites if s != nearest_neighbor)
        return determine_starting_bisector(w, next_r, width, intersection['point'], find_bisector)
    
    else:
        w = find_correct_w(w, nearest_neighbor, find_bisector)
        starting_bisector = find_bisector(w, nearest_neighbor)
        
        return {
            'starting_bisector': starting_bisector,
            'w': w,
            'nearest_neighbor': nearest_neighbor,
            'starting_intersection': intersection['point'] if intersection else w.site
        }


def find_correct_w(w: Site, nearest_neighbor: Site, find_bisector):
    """Ensure that the starting point is correct and would not result in a trapped bisector"""
    starting_bisector = find_bisector(w, nearest_neighbor)
    
    w_trap = None
    min_dist = float('inf')
    
    for bisector in w.bisectors:
        hop_to = find_hop_to(bisector, w)
        if is_bisector_trapped(hop_to, starting_bisector):
            dist = manhattan_distance(hop_to.site, nearest_neighbor.site)
            if dist < min_dist:
                min_dist = dist
                w_trap = hop_to
    
    if w_trap:
        return find_correct_w(w_trap, nearest_neighbor, find_bisector)
    else:
        return w


def check_for_orphans(trapper: Site, trapped: Site, go_up: bool, find_bisector):
    """Recursively check for orphaned bisectors"""
    orphans = []
    
    for bisector in trapped.bisectors:
        hop_to = find_hop_to(bisector, trapped)
        if (go_up == (hop_to.site[1] < trapped.site[1]) and 
            is_bisector_trapped(trapper, bisector)):
            
            merge_line = find_bisector(hop_to, trapper)
            extreme = get_extreme_point(merge_line, go_up)
            orphans.append({'bisector': bisector, 'hop_to': hop_to, 'extreme': extreme})
    
    if not orphans:
        return None
    
    # Sort by extreme point
    orphans.sort(key=lambda x: x['extreme'], reverse=go_up)
    return orphans[0]['bisector']


def curry_find_bisector(callback, width: int, height: int):
    """Curry find bisector function with the current width, height"""
    return lambda p1, p2: callback(p1, p2, width, height)


def find_l1_bisector(p1: Site, p2: Site, width: int, height: int):
    """Generate L1 bisector between two sites"""
    x_distance = p1.site[0] - p2.site[0]
    y_distance = p1.site[1] - p2.site[1]
    
    midpoint = [
        (p1.site[0] + p2.site[0]) / 2,
        (p1.site[1] + p2.site[1]) / 2
    ]
    
    # Handle special cases first
    if abs(x_distance) == abs(y_distance):
        raise ValueError(f"Square bisector: Points {p1.site} and {p2.site} are points on a square")
    
    if same_point(p1.site, p2.site):
        raise ValueError(f"Duplicate point: Points {p1.site} and {p2.site} are duplicates")
    
    # Vertical bisector (horizontal line)
    if abs(x_distance) == 0:
        points = [
            [0, midpoint[1]],
            [width, midpoint[1]]
        ]
        return Bisector([p1, p2], False, points, [], False)
    
    # Horizontal bisector (vertical line)
    if abs(y_distance) == 0:
        points = [
            [midpoint[0], 0],
            [midpoint[0], height]
        ]
        return Bisector([p1, p2], True, points, [], False)
    
    # Diagonal bisector
    slope = -1 if y_distance / x_distance > 0 else 1
    intercept = midpoint[1] - midpoint[0] * slope
    
    if abs(x_distance) >= abs(y_distance):
        # More horizontal than vertical - bisector goes up/down
        vertices = [
            [max(0, min(width, (p1.site[1] - intercept) / slope)), p1.site[1]],
            [max(0, min(width, (p2.site[1] - intercept) / slope)), p2.site[1]]
        ]
        up = True
        
        sorted_verts = sorted(vertices, key=lambda v: v[1])
        points = sorted(
            [[sorted_verts[0][0], 0]] + sorted_verts + [[sorted_verts[1][0], height]],
            key=lambda p: p[1]
        )
    else:
        # More vertical than horizontal - bisector goes left/right
        vertices = [
            [p1.site[0], max(0, min(height, (p1.site[0] * slope) + intercept))],
            [p2.site[0], max(0, min(height, (p2.site[0] * slope) + intercept))]
        ]
        up = False
        
        sorted_verts = sorted(vertices, key=lambda v: v[0])
        points = sorted(
            [[0, sorted_verts[0][1]]] + sorted_verts + [[width, sorted_verts[1][1]]],
            key=lambda p: p[0]
        )
    
    return Bisector([p1, p2], up, points, [], False)


def clear_out_orphans(orphanage: Site, trap_point: Site):
    """Clear out orphans when a new merge line is created"""
    return [b for b in orphanage.bisectors if not is_bisector_trapped(trap_point, b)]


def find_hop_to(bisector: Bisector, hop_from: Site) -> Site:
    """Find other point across a bisector"""
    return next(s for s in bisector.sites if s != hop_from)


def manhattan_distance(p1: List[float], p2: List[float]) -> float:
    """Find L1 distance"""
    return abs(p1[0] - p2[0]) + abs(p1[1] - p2[1])


def is_bisector_trapped(trap_point: Site, bisector: Bisector) -> bool:
    """Determine if bisector is trapped in a site's polygon"""
    return all(manhattan_distance(trap_point.site, point) <= 
              min(manhattan_distance(bisector.sites[0].site, point),
                  manhattan_distance(bisector.sites[1].site, point))
              for point in bisector.points)


def get_extreme_point(bisector: Bisector, go_up: bool) -> float:
    """Find the highest or lowest point of a potential bisector"""
    if go_up:
        return max(p[1] for p in bisector.points)
    else:
        return min(p[1] for p in bisector.points)


def trim_bisector(target: Bisector, intersector: Optional[Bisector], 
                 intersection: List[float]):
    """Trim a bisector at a particular point"""
    if not intersector:
        return
    
    polygon_site = next(s for s in intersector.sites if s not in target.sites)
    
    new_points = [p for p in target.points 
                  if (manhattan_distance(p, target.sites[0].site) < 
                      manhattan_distance(p, polygon_site.site) and
                      manhattan_distance(p, target.sites[1].site) < 
                      manhattan_distance(p, polygon_site.site))]
    
    new_points.append(intersection)
    
    target.points = sorted(new_points, key=lambda p: p[1] if target.up else p[0])


def is_new_bisector_upward(hop_to: Site, hop_from: Site, site: Site, go_up: bool) -> bool:
    """Check to see if a bisector is traveling upward or downward"""
    if hop_to.site[0] == site.site[0]:
        return site.site[1] > hop_to.site[1]
    
    slope = (hop_to.site[1] - site.site[1]) / (hop_to.site[0] - site.site[0])
    intercept = hop_to.site[1] - (slope * hop_to.site[0])
    
    is_above_line = hop_from.site[1] > (slope * hop_from.site[0]) + intercept
    return is_above_line


def bisector_intersection(b1: Bisector, b2: Bisector) -> Optional[List[float]]:
    """Find intersection of two bisectors, if it exists"""
    if b1 == b2 or not b1.points or not b2.points:
        return None
    
    # Ensure we have at least 2 points in each bisector
    if len(b1.points) < 2 or len(b2.points) < 2:
        return None
    
    for i in range(len(b1.points) - 1):
        for j in range(len(b2.points) - 1):
            intersect = segment_intersection(
                [b1.points[i], b1.points[i+1]], 
                [b2.points[j], b2.points[j+1]]
            )
            if intersect:
                return intersect
    
    return None


def segment_intersection(l1: List[List[float]], l2: List[List[float]]) -> Optional[List[float]]:
    """Find intersection of two line segments"""
    denom = ((l2[1][1] - l2[0][1]) * (l1[1][0] - l1[0][0]) - 
             (l2[1][0] - l2[0][0]) * (l1[1][1] - l1[0][1]))
    
    if denom == 0:
        return None
    
    ua = ((l2[1][0] - l2[0][0]) * (l1[0][1] - l2[0][1]) - 
          (l2[1][1] - l2[0][1]) * (l1[0][0] - l2[0][0])) / denom
    
    ub = ((l1[1][0] - l1[0][0]) * (l1[0][1] - l2[0][1]) - 
          (l1[1][1] - l1[0][1]) * (l1[0][0] - l2[0][0])) / denom
    
    if not (0 <= ua <= 1 and 0 <= ub <= 1):
        return None
    
    return [
        l1[0][0] + ua * (l1[1][0] - l1[0][0]),
        l1[0][1] + ua * (l1[1][1] - l1[0][1])
    ]


def same_point(p1: List[float], p2: List[float]) -> bool:
    """Determine if two points are the same point"""
    return p1[0] == p2[0] and p1[1] == p2[1]


def visualize_voronoi(sites: List[Site], width: int, height: int, 
                     merge_process: bool = False):
    """Visualize the Voronoi diagram using matplotlib"""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
    
    # Main diagram
    ax1.set_title("Voronoi Diagram")
    ax1.set_xlim(0, width)
    ax1.set_ylim(0, height)
    ax1.set_aspect('equal')
    
    for site in sites:
        # Draw polygon
        if site.polygon_points:
            polygon = patches.Polygon(site.polygon_points, closed=True, 
                                     fill=True, facecolor='lightblue', 
                                     edgecolor='black', linewidth=1)
            ax1.add_patch(polygon)
        
        # Draw site point
        ax1.plot(site.site[0], site.site[1], 'ko', markersize=3)
    
    # Merge process
    if merge_process:
        ax2.set_title("Merge Process")
        ax2.set_xlim(0, width)
        ax2.set_ylim(0, height)
        ax2.set_aspect('equal')
        
        colors = {
            4: '#4286f4',
            8: '#44f453',
            16: '#931d78',
            32: '#ff3c35',
            64: '#f4ad42',
            128: '#009182',
            256: '#993300',
            512: '#669999',
            1024: '#800000',
            2048: '#333300'
        }
        
        for site in sites:
            # Draw bisectors
            for bisector in site.bisectors:
                if bisector.points:
                    points = np.array(bisector.points)
                    color = colors.get(bisector.merge_line, '#000000') if bisector.merge_line else '#000000'
                    ax2.plot(points[:, 0], points[:, 1], color=color, linewidth=1)
            
            # Draw site point
            ax2.plot(site.site[0], site.site[1], 'ko', markersize=3)
    
    plt.tight_layout()
    plt.show()


def main():
    """Main function to generate and visualize L1 Voronoi diagram"""
    width = 400
    height = 400
    num_points = 8  # Use fewer points for testing
    
    # Generate random points with normal distribution
    def random_normal(sharpness):
        return sum(random.random() for _ in range(sharpness)) / sharpness
    
    raw_points = [
        [int(random_normal(2) * width), int(random_normal(2) * height)]
        for _ in range(num_points)
    ]
    
    # Sort points for consistency
    sites = sorted(raw_points, key=lambda p: (p[0], p[1]))
    
    print("Generated points:")
    for point in sites:
        print(f"  {point}")
    
    try:
        # Generate Voronoi diagram
        voronoi_sites = generate_l1_voronoi(raw_points, width, height, True)
        
        # Visualize
        visualize_voronoi(voronoi_sites, width, height, merge_process=True)
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()