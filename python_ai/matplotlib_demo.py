"""
Matplotlib demo for L1 (Manhattan) Voronoi diagram
"""

from voronoi import generate_l1_voronoi
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import random

def random_normal(sharpness):
    """Generate a normally distributed random number"""
    return sum(random.random() for _ in range(sharpness)) / sharpness

def create_custom_voronoi_demo():
    """Create a custom Voronoi visualization with matplotlib"""
    
    # Set up the figure
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    fig.suptitle('L1 (Manhattan) Voronoi Diagram Demonstration', fontsize=16)
    
    # Demo configurations
    configs = [
        {'points': 4, 'title': 'Simple 4-Point Configuration'},
        {'points': 8, 'title': '8-Point Configuration'},
        {'points': 16, 'title': '16-Point Configuration'},
        {'points': 32, 'title': '32-Point Configuration'}
    ]
    
    for idx, (ax, config) in enumerate(zip(axes.flat, configs)):
        width, height = 200, 200
        num_points = config['points']
        
        # Generate random points
        raw_points = [
            [int(random_normal(2) * width), int(random_normal(2) * height)]
            for _ in range(num_points)
        ]
        
        # Generate Voronoi diagram
        voronoi_sites = generate_l1_voronoi(raw_points, width, height, nudge_data=True)
        
        # Set up the subplot
        ax.set_title(config['title'])
        ax.set_xlim(0, width)
        ax.set_ylim(0, height)
        ax.set_aspect('equal')
        ax.grid(True, alpha=0.3)
        
        # Create a colormap for the cells
        colors = plt.cm.Set3(np.linspace(0, 1, len(voronoi_sites)))
        
        # Draw Voronoi cells
        for i, site in enumerate(voronoi_sites):
            if site.polygon_points:
                # Draw filled polygon
                polygon = patches.Polygon(
                    site.polygon_points, 
                    closed=True, 
                    fill=True, 
                    facecolor=colors[i], 
                    edgecolor='black', 
                    linewidth=1.5,
                    alpha=0.7
                )
                ax.add_patch(polygon)
                
                # Draw site point
                ax.plot(site.site[0], site.site[1], 'ko', markersize=6)
                
                # Add site label
                ax.text(site.site[0], site.site[1], f'{i+1}', 
                       ha='center', va='center', fontsize=8, 
                       color='white', weight='bold')
        
        # Add statistics text
        ax.text(0.02, 0.98, f'Points: {num_points}\nCells: {len(voronoi_sites)}', 
               transform=ax.transAxes, fontsize=9,
               verticalalignment='top',
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig('matplotlib_voronoi_demo.png', dpi=150, bbox_inches='tight')
    print("Matplotlib demo saved to 'matplotlib_voronoi_demo.png'")
    
    # Create a single detailed view
    fig2, ax = plt.subplots(1, 1, figsize=(10, 10))
    
    # Use more points for the detailed view
    width, height = 400, 400
    raw_points = [
        [int(random_normal(2) * width), int(random_normal(2) * height)]
        for _ in range(25)
    ]
    
    voronoi_sites = generate_l1_voronoi(raw_points, width, height, nudge_data=True)
    
    ax.set_title('Detailed L1 Voronoi Diagram (25 points)', fontsize=14)
    ax.set_xlim(0, width)
    ax.set_ylim(0, height)
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.2)
    
    # Use a different colormap
    colors = plt.cm.tab20(np.linspace(0, 1, len(voronoi_sites)))
    
    # Draw cells with more detail
    for i, site in enumerate(voronoi_sites):
        if site.polygon_points:
            # Draw filled polygon
            polygon = patches.Polygon(
                site.polygon_points, 
                closed=True, 
                fill=True, 
                facecolor=colors[i % 20], 
                edgecolor='black', 
                linewidth=1,
                alpha=0.6
            )
            ax.add_patch(polygon)
            
            # Draw site point with larger marker
            ax.plot(site.site[0], site.site[1], 'o', 
                   color='black', markersize=8, 
                   markeredgecolor='white', markeredgewidth=1)
    
    # Add legend with point coordinates
    legend_elements = []
    for i, site in enumerate(voronoi_sites[:10]):  # Show first 10 in legend
        legend_elements.append(
            plt.Line2D([0], [0], marker='o', color='w', 
                      markerfacecolor=colors[i % 20], 
                      markersize=10, 
                      label=f'Site {i+1}: ({site.site[0]}, {site.site[1]})')
        )
    
    if len(voronoi_sites) > 10:
        legend_elements.append(
            plt.Line2D([0], [0], marker='o', color='w', 
                      markerfacecolor='gray', 
                      markersize=10, 
                      label=f'... and {len(voronoi_sites) - 10} more')
        )
    
    ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1, 1))
    
    plt.tight_layout()
    plt.savefig('matplotlib_voronoi_detailed.png', dpi=150, bbox_inches='tight')
    print("Detailed visualization saved to 'matplotlib_voronoi_detailed.png'")
    
    # Print statistics
    print(f"\nDetailed View Statistics:")
    print(f"  Total points: {len(raw_points)}")
    print(f"  Voronoi cells: {len(voronoi_sites)}")
    print(f"  Average vertices per cell: {np.mean([len(s.polygon_points) for s in voronoi_sites]):.1f}")
    print(f"  Average bisectors per cell: {np.mean([len(s.bisectors) for s in voronoi_sites]):.1f}")

if __name__ == "__main__":
    try:
        print("Starting matplotlib Voronoi demo...")
        create_custom_voronoi_demo()
        print("Demo completed successfully!")
    except Exception as e:
        print(f"Error running demo: {e}")
        import traceback
        traceback.print_exc()