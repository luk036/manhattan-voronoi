def show_bisectors_only():
    """Display only the bisector lines"""
    
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
    
    width, height = 400, 400
    
    print("Generating L1 Voronoi diagram...")
    print("Points:")
    for i, p in enumerate(points):
        print(f"  {i+1}. [{p[0]:3d}, {p[1]:3d}]")
    
    # Generate Voronoi diagram
    sites = generate_l1_voronoi(points, width, height, nudge_data=True)
    
    print(f"\nGenerated {len(sites)} Voronoi cells")
    
    # Create figure
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    ax.set_title('L1 Voronoi Bisectors\n(only 0°, 45°, 90°, and 180° lines)', fontsize=14)
    ax.set_xlim(0, width)
    ax.set_ylim(0, height)
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)
    
    # Collect all unique bisectors
    all_bisectors = []
    seen = set()
    for site in sites:
        for bisector in site.bisectors:
            # Create unique key from site coordinates
            key = tuple(sorted([tuple(bisector.sites[0].site), tuple(bisector.sites[1].site)]))
            if key not in seen:
                seen.add(key)
                all_bisectors.append(bisector)
    
    # Draw each bisector
    colors = plt.cm.rainbow(np.linspace(0, 1, len(all_bisectors)))
    
    line_stats = []
    for i, bisector in enumerate(all_bisectors):
        if bisector.points and len(bisector.points) >= 2:
            points = np.array(bisector.points)
            ax.plot(points[:, 0], points[:, 1], color=colors[i], 
                    linewidth=3, alpha=0.8, label=f'Bisector {i+1}')
            
            # Calculate line angle
            p1 = bisector.points[0]
            p2 = bisector.points[1] if len(bisector.points) > 1 else bisector.points[-1]
            dx = p2[0] - p1[0]
            dy = p2[1] - p1[1]
            
            if dx != 0 or dy != 0:
                angle_deg = math.degrees(math.atan2(dy, dx)) % 180
                line_stats.append(angle_deg)
                
                # Determine line type
                if abs(angle_deg - 0) < 1 or abs(angle_deg - 180) < 1:
                    line_type = "Horizontal"
                elif abs(angle_deg - 90) < 1:
                    line_type = "Vertical"
                elif abs(angle_deg - 45) < 1 or abs(angle_deg - 135) < 1:
                    line_type = "Diagonal"
                else:
                    line_type = f"Other ({angle_deg:.1f}°)"
                
                # Add label at midpoint
                mid_x = (p1[0] + p2[0]) / 2
                mid_y = (p1[1] + p2[1]) / 2
                ax.text(mid_x, mid_y, f'{angle_deg:.0f}°\n{line_type}', 
                        fontsize=9, ha='center', va='center',
                        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))
    
    # Draw site points
    for i, site in enumerate(sites):
        ax.plot(site.site[0], site.site[1], 'ko', markersize=10)
        ax.text(site.site[0], site.site[1] + 15, f'S{i+1}\n({site.site[0]}, {site.site[1]})', 
                ha='center', va='bottom', fontsize=10, weight='bold')
    
    # Add legend
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    plt.tight_layout()
    plt.savefig('bisectors_only.png', dpi=150, bbox_inches='tight')
    print("\nVisualization saved to 'bisectors_only.png'")
    
    # Print statistics
    print("\nBisector Statistics:")
    print(f"  Total bisectors: {len(all_bisectors)}")
    
    # Count lines by type
    horizontal = sum(1 for a in line_stats if abs(a) < 1 or abs(a - 180) < 1)
    vertical = sum(1 for a in line_stats if abs(a - 90) < 1)
    diagonal = sum(1 for a in line_stats if abs(a - 45) < 1 or abs(a - 135) < 1)
    
    print(f"  Horizontal lines (0°/180°): {horizontal}")
    print(f"  Vertical lines (90°): {vertical}")
    print(f"  Diagonal lines (45°/135°): {diagonal}")
    
    if horizontal + vertical + diagonal < len(line_stats):
        other = len(line_stats) - horizontal - vertical - diagonal
        print(f"  Other angles: {other}")
        print(f"  Angle values: {[round(a, 1) for a in line_stats if not (abs(a) < 1 or abs(a - 90) < 1 or abs(a - 45) < 1 or abs(a - 135) < 1)]}")
    
    # Verify all lines are at expected angles
    expected_angles = [0, 45, 90, 135, 180]
    unexpected = [a for a in line_stats if all(abs(a - e) > 2 for e in expected_angles)]
    if unexpected:
        print(f"\nWARNING: Found unexpected angles: {[round(a, 1) for a in unexpected]}")
    else:
        print("\nSUCCESS: All bisectors are at expected angles (0°, 45°, 90°, 135°, or 180°)")

if __name__ == "__main__":
    import matplotlib
    from voronoi_clean import generate_l1_voronoi
    import matplotlib.pyplot as plt
    import numpy as np
    import math
    matplotlib.use('Agg')  # Use non-interactive backend
    show_bisectors_only()