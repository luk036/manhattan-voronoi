"""
Demo script for L1 (Manhattan) Voronoi diagram generation
"""

from voronoi import generate_l1_voronoi, visualize_voronoi
import random
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt

def random_normal(sharpness):
    """Generate a normally distributed random number"""
    return sum(random.random() for _ in range(sharpness)) / sharpness

def generate_random_points(num_points, width, height):
    """Generate random points with normal distribution"""
    return [
        [int(random_normal(2) * width), int(random_normal(2) * height)]
        for _ in range(num_points)
    ]

def main():
    """Main demo function"""
    width = 400
    height = 400
    num_points = 16  # Start with fewer points for stability
    
    print(f"Generating L1 Voronoi diagram with {num_points} points...")
    
    # Generate random points
    raw_points = generate_random_points(num_points, width, height)
    
    # Sort points for consistency (as required by the algorithm)
    sites = sorted(raw_points, key=lambda p: (p[0], p[1]))
    
    print("\nGenerated points (sorted by x, then y):")
    for i, point in enumerate(sites):
        print(f"  {i+1}. [{point[0]:3d}, {point[1]:3d}]")
    
    try:
        # Generate Voronoi diagram with nudge data to avoid square bisectors
        voronoi_sites = generate_l1_voronoi(raw_points, width, height, nudge_data=True)
        
        print(f"\nSuccessfully generated Voronoi diagram with {len(voronoi_sites)} cells")
        
        # Create visualization
        plt.figure(figsize=(12, 6))
        visualize_voronoi(voronoi_sites, width, height, merge_process=True)
        plt.savefig('manhattan_voronoi.png', dpi=150, bbox_inches='tight')
        print("\nVisualization saved to 'manhattan_voronoi.png'")
        
        # Show some statistics
        print("\nVoronoi cell statistics:")
        for i, site in enumerate(voronoi_sites[:5]):  # Show first 5 cells
            print(f"  Cell {i+1}: Site at {site.site}, "
                  f"{len(site.bisectors)} bisectors, "
                  f"{len(site.polygon_points)} polygon vertices")
        
        if len(voronoi_sites) > 5:
            print(f"  ... and {len(voronoi_sites) - 5} more cells")
            
    except Exception as e:
        print(f"\nError generating Voronoi diagram: {e}")
        print("\nThis might happen with certain point configurations.")
        print("Try running the script again with different random points.")
        
        # Try with a simpler configuration as fallback
        print("\nTrying with a simple 4-point configuration...")
        simple_points = [
            [50, 50],
            [150, 80],
            [120, 150],
            [80, 120]
        ]
        
        voronoi_sites = generate_l1_voronoi(simple_points, width, height, nudge_data=True)
        visualize_voronoi(voronoi_sites, width, height, merge_process=True)
        plt.savefig('manhattan_voronoi_simple.png', dpi=150, bbox_inches='tight')
        print("Simple visualization saved to 'manhattan_voronoi_simple.png'")

if __name__ == "__main__":
    main()