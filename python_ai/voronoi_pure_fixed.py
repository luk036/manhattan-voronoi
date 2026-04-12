 """
Pure L1 Voronoi implementation that only draws the actual bisector lines
"""

from voronoi_clean import generate_l1_voronoi
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import math

def visualize_pure_voronoi(sites, width, height):
    """Visualize only the actual bisector lines without polygon filling"""
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
    
    # Left plot: Show only bisector lines
    ax1.set_title('Pure L1 Voronoi Bisectors\n(only 0°, 45°, 90°, 180° lines)', fontsize=14)
    ax1.set_xlim(0, width)
    ax1.set_ylim(0, height)
    ax1.set_aspect('equal')
    ax1.grid(True, alpha=0.3)
    
    # Collect all unique bisectors
    all_bisectors = set()
    for site in sites:
        for bisector in site.bisectors:
            # Create a unique identifier for each bisector (pair of sites)
            site_ids = sorted([id(bisector.sites[0]), id(bisector.sites[1])])
            all_bisectors.add((site_ids[0], site_ids[1], bisector))
    
    # Draw each bisector line
    colors = plt.cm.tab10(np.linspace(0, 1, len(all_bisectors)))
    for i, (_, _, bisector) in enumerate(all_bisectors):
        if bisector.points and len(bisector.points) >= 2:
            points = np.array(bisector.points)
            ax1.plot(points[:, 0], points[:, 1], color=colors[i % 10], 
                    linewidth=2, alpha=0.8)
            
            # Calculate and display line angle
            p1 = bisector.points[0]
            p2 = bisector.points[1]
            dx = p2[0] - p1[0]
            dy = p2[1] - p1[1]
            if dx != 0 or dy != 0:
                angle_deg = math.degrees(math.atan2(dy, dx))
                # Normalize to 0-180 degrees
                angle_deg = angle_deg % 180
                
                # Determine line type
                if abs(angle_deg - 0) < 1 or abs(angle_deg - 180) < 1:
                    line_type = "Horizontal (0°/180°)"
                elif abs(angle_deg - 90) < 1:
                    line_type = "Vertical (90°)"
                elif abs(angle_deg - 45) < 1 or abs(angle_deg - 135) < 1:
                    line_type = "Diagonal (45°/135°)"
                else:
                    line_type = f"Other ({angle_deg:.1f}°)"
                
                # Add label at midpoint
                mid_x = (p1[0] + p2[0]) / 2
                mid_y = (p1[1] + p2[1]) / 2
                ax1.text(mid_x, mid_y, f'{angle_deg:.0f}°', 
                        fontsize=8, ha='center', va='bottom',
                        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.7))
    
    # Draw site points
    for i, site in enumerate(sites):
        ax1.plot(site.site[0], site.site[1], 'ko', markersize=8)
        ax1.text(site.site[0], site.site[1] + 10, f'S{i+1}', 
                ha='center', va='bottom', fontsize=10, weight='bold')
    
    # Right plot: Show Voronoi cells with proper construction
    ax2.set_title('L1 Voronoi Cells\n(constructed from bisector intersections)', fontsize=14)
    ax2.set_xlim(0, width)
    ax2.set_ylim(0, height)
    ax2.set_aspect('equal')
    ax2.grid(True, alpha=0.3)
    
    # Construct Voronoi cells properly
    colors = plt.cm.Set3(np.linspace(0, 1, len(sites)))
    
    for i, site in enumerate(sites):
        # Find all points that belong to this site's Voronoi cell
        cell_points = []
        
        # Add all bisector points that are closer to this site
        for bisector in site.bisectors:
            for point in bisector.points:
                # Check if this point is on the boundary of this site's cell
                other_site = bisector.sites[0] if bisector.sites[1] == site else bisector.sites[1]
                if abs(manhattan_distance(point, site.site) - manhattan_distance(point, other_site.site)) < 0.1:
                    # Check if no other site is closer
                    is_boundary = True
                    for check_site in sites:
                        if check_site != site and check_site != other_site:
                            if manhattan_distance(point, check_site.site) < manhattan_distance(point, site.site):
                                is_boundary = False
                                break
                    if is_boundary and point not in cell_points:
                        cell_points.append(point)
        
        # Add intersections between this site's bisectors
        for j, b1 in enumerate(site.bisectors):
            for b2 in site.bisectors[j+1:]:
                intersection = bisector_intersection(b1, b2)
                if intersection:
                    # Check if this intersection point is closest to our site
                    is_closest = True
                    for other_site in sites:
                        if other_site != site:
                            if manhattan_distance(intersection, other_site.site) < manhattan_distance(intersection, site.site):
                                is_closest = False
                                break
                    if is_closest and intersection not in cell_points:
                        cell_points.append(intersection)
        
        # Sort points by angle around the site
        if cell_points:
            cell_points.sort(key=lambda p: angle(site.site, p))
            
            # Draw filled polygon
            polygon = patches.Polygon(
                cell_points,
                closed=True,
                fill=True,
                facecolor=colors[i],
                edgecolor='black',
                linewidth=2,
                alpha=0.6
            )
            ax2.add_patch(polygon)
        
        # Draw site point
        ax2.plot(site.site[0], site.site[1], 'ko', markersize=8)
        ax2.text(site.site[0], site.site[1], f'{i+1}', 
                ha='center', va='center', fontsize=10, 
                color='white', weight='bold')
    
    plt.tight_layout()
    plt.show()
    
    # Print line statistics
    print("\nBisector Line Statistics:")
    line_angles = []
    for _, _, bisector in all_bisectors:
        if bisector.points and len(bisector.points) >= 2:
            p1 = bisector.points[0]
            p2 = bisector.points[1]
            dx = p2[0] - p1[0]
            dy = p2[1] - p1[1]
            if dx != 0 or dy != 0:
                angle_deg = math.degrees(math.atan2(dy, dx)) % 180
                line_angles.append(round(angle_deg, 1))
    
    unique_angles = sorted(set(line_angles))
    print(f"  Total bisectors: {len(all_bisectors)}")
    print(f"  Unique angles: {unique_angles}")
    
    # Count lines by type
    horizontal = sum(1 for a in line_angles if abs(a) < 1 or abs(a - 180) < 1)
    vertical = sum(1 for a in line_angles if abs(a - 90) < 1)
    diagonal_45 = sum(1 for a in line_angles if abs(a - 45) < 1 or abs(a - 135) < 1)
    
    print(f"  Horizontal lines (0°/180°): {horizontal}")
    print(f"  Vertical lines (90°): {vertical}")
    print(f"  Diagonal lines (45°/135°): {diagonal_45}")
    if horizontal + vertical + diagonal_45 < len(line_angles):
        print(f"  Other angles: {len(line_angles) - horizontal - vertical - diagonal_45}")


def manhattan_distance(p1, p2):
    """Calculate Manhattan distance"""
    return abs(p1[0] - p2[0]) + abs(p1[1] - p2[1])


def angle(p1, p2):
    """Calculate angle between two points"""
    angle = math.atan2(p2[1] - p1[1], p2[0] - p1[0])
    if angle < 0:
        angle = math.pi + math.pi + angle
    return angle


def bisector_intersection(b1, b2):
    """Find intersection of two bisectors"""
    if b1 == b2 or not b1.points or not b2.points:
        return None
    
    if len(b1.points) < 2 or len(b2.points) < 2:
        return None
    
    for i in range(len(b1.points) - 1):
        for j in range(len(b2.points) - 1):
            # Simple line intersection
            x1, y1 = b1.points[i]
            x2, y2 = b1.points[i+1]
            x3, y3 = b2.points[j]
            x4, y4 = b2.points[j+1]
            
            denom = (x1-x2)*(y3-y4) - (y1-y2)*(x3-x4)
            if abs(denom) > 1e-10:
                t = ((x1-x3)*(y3-y4) - (y1-y3)*(x3-x4)) / denom
                u = -((x1-x2)*(y1-y3) - (y1-y2)*(x1-x3)) / denom
                
                if 0 <= t <= 1 and 0 <= u <= 1:
                    ix = x1 + t * (x2 - x1)
                    iy = y1 + t * (y2 - y1)
                    return [ix, iy]
    
    return None


def main():
    """Main function"""
    width, height = 400, 400
    
    # Use points that won't create square bisectors
    points = [
        [50, 50],
        [150, 80],
        [120, 150],
        [80, 120],
        [200, 100],
        [250, 200],
        [300, 250],
        [100, 250]
    ]
    
    print("Generating pure L1 Voronoi diagram...")
    print("Points:")
    for i, p in enumerate(points):
        print(f"  {i+1}. [{p[0]:3d}, {p[1]:3d}]")
    
    # Generate Voronoi diagram
    sites = generate_l1_voronoi(points, width, height, nudge_data=True)
    
    print(f"\nGenerated {len(sites)} Voronoi cells")
    
    # Visualize with pure bisector lines
    import matplotlib
    matplotlib.use('TkAgg')  # Use interactive backend
    visualize_pure_voronoi(sites, width, height)


if __name__ == "__main__":
    main()
