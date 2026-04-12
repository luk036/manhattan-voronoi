"""
Compare original and clean L1 Voronoi implementations
"""

from voronoi import generate_l1_voronoi as original_voronoi
from voronoi_clean import generate_l1_voronoi as clean_voronoi
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

def compare_implementations():
    """Compare original and clean Voronoi implementations"""
    
    # Set up the figure
    fig, axes = plt.subplots(1, 2, figsize=(16, 8))
    fig.suptitle('L1 Voronoi Diagram Comparison', fontsize=16)
    
    # Use the same points for both implementations
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
    
    width, height = 400, 400
    
    # Generate Voronoi diagrams
    try:
        original_sites = original_voronoi(points, width, height, nudge_data=True)
    except ValueError as e:
        print(f"Original implementation failed: {e}")
        print("Using modified points for original implementation...")
        # Slightly modify points to avoid square bisectors
        modified_points = []
        for i, p in enumerate(points):
            modified_points.append([p[0] + 0.1 if i % 2 == 0 else p[0], 
                                  p[1] + 0.1 if i % 2 == 1 else p[1]])
        original_sites = original_voronoi(modified_points, width, height, nudge_data=True)
    
    clean_sites = clean_voronoi(points, width, height, nudge_data=True)
    
    # Plot original implementation
    ax1 = axes[0]
    ax1.set_title('Original Implementation\n(with extra lines)', fontsize=12)
    ax1.set_xlim(0, width)
    ax1.set_ylim(0, height)
    ax1.set_aspect('equal')
    ax1.grid(True, alpha=0.3)
    
    colors1 = plt.cm.Set3(np.linspace(0, 1, len(original_sites)))
    
    for i, site in enumerate(original_sites):
        if site.polygon_points:
            polygon = patches.Polygon(
                site.polygon_points, 
                closed=True, 
                fill=True, 
                facecolor=colors1[i], 
                edgecolor='black', 
                linewidth=1,
                alpha=0.7
            )
            ax1.add_patch(polygon)
            ax1.plot(site.site[0], site.site[1], 'ko', markersize=6)
    
    # Plot clean implementation
    ax2 = axes[1]
    ax2.set_title('Clean Implementation\n(only 0°, 45°, 90°, 180° lines)', fontsize=12)
    ax2.set_xlim(0, width)
    ax2.set_ylim(0, height)
    ax2.set_aspect('equal')
    ax2.grid(True, alpha=0.3)
    
    colors2 = plt.cm.Set3(np.linspace(0, 1, len(clean_sites)))
    
    for i, site in enumerate(clean_sites):
        if site.polygon_points:
            polygon = patches.Polygon(
                site.polygon_points, 
                closed=True, 
                fill=True, 
                facecolor=colors2[i], 
                edgecolor='black', 
                linewidth=2,
                alpha=0.7
            )
            ax2.add_patch(polygon)
            ax2.plot(site.site[0], site.site[1], 'ko', markersize=6)
            ax2.text(site.site[0], site.site[1], f'{i+1}', 
                    ha='center', va='center', fontsize=8, 
                    color='white', weight='bold')
    
    # Add statistics
    stats_text1 = f"Original:\n"
    stats_text1 += f"  Cells: {len(original_sites)}\n"
    stats_text1 += f"  Avg vertices: {np.mean([len(s.polygon_points) for s in original_sites]):.1f}\n"
    stats_text1 += f"  Max vertices: {max([len(s.polygon_points) for s in original_sites])}"
    
    stats_text2 = f"Clean:\n"
    stats_text2 += f"  Cells: {len(clean_sites)}\n"
    stats_text2 += f"  Avg vertices: {np.mean([len(s.polygon_points) for s in clean_sites]):.1f}\n"
    stats_text2 += f"  Max vertices: {max([len(s.polygon_points) for s in clean_sites])}"
    
    ax1.text(0.02, 0.98, stats_text1, transform=ax1.transAxes, fontsize=10,
             verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    ax2.text(0.02, 0.98, stats_text2, transform=ax2.transAxes, fontsize=10,
             verticalalignment='top', bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig('voronoi_comparison.png', dpi=150, bbox_inches='tight')
    print("Comparison saved to 'voronoi_comparison.png'")
    
    # Print detailed statistics
    print("\nDetailed Statistics:")
    print("\nOriginal Implementation:")
    for i, site in enumerate(original_sites):
        print(f"  Cell {i+1}: {len(site.polygon_points)} vertices")
    
    print("\nClean Implementation:")
    for i, site in enumerate(clean_sites):
        print(f"  Cell {i+1}: {len(site.polygon_points)} vertices")
    
    # Check line angles in clean implementation
    print("\nLine angles in clean implementation:")
    for i, site in enumerate(clean_sites):
        if site.polygon_points and len(site.polygon_points) > 1:
            angles = []
            for j in range(len(site.polygon_points)):
                p1 = site.polygon_points[j]
                p2 = site.polygon_points[(j+1) % len(site.polygon_points)]
                dx = p2[0] - p1[0]
                dy = p2[1] - p1[1]
                if dx != 0 or dy != 0:
                    angle_deg = math.degrees(math.atan2(dy, dx))
                    # Normalize to 0-180 degrees
                    angle_deg = angle_deg % 180
                    angles.append(round(angle_deg, 1))
            unique_angles = sorted(set(angles))
            print(f"  Cell {i+1}: {unique_angles}")

if __name__ == "__main__":
    import math
    compare_implementations()
