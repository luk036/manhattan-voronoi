from voronoi import generate_l1_voronoi, visualize_voronoi

def test_simple():
    """Test with a simple 4-point configuration"""
    width = 100
    height = 100
    
    # Use 4 simple points avoiding square bisectors
    raw_points = [
        [20, 25],  # Bottom-left
        [75, 30],  # Bottom-right
        [70, 75],  # Top-right
        [25, 70],  # Top-left
    ]
    
    print("Test points:")
    for point in raw_points:
        print(f"  {point}")
    
    try:
        # Generate Voronoi diagram
        voronoi_sites = generate_l1_voronoi(raw_points, width, height, True)
        print(f"Generated {len(voronoi_sites)} Voronoi sites")
        
        # Debug: Check polygon points
        for i, site in enumerate(voronoi_sites):
            print(f"Site {i}: {site.site}")
            print(f"  Polygon points: {site.polygon_points}")
            print(f"  Bisectors: {len(site.bisectors)}")
        
        # Visualize and save
        import matplotlib
        matplotlib.use('Agg')  # Use non-interactive backend
        visualize_voronoi(voronoi_sites, width, height, merge_process=True)
        import matplotlib.pyplot as plt
        plt.savefig('voronoi_test.png', dpi=150, bbox_inches='tight')
        print("Saved visualization to voronoi_test.png")
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    test_simple()