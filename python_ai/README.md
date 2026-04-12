# Manhattan (L1) Voronoi Diagram in Python

This is a Python implementation of the rectilinear (L1) Voronoi diagram using Lee and Wong's divide-and-conquer algorithm. The code is converted from the original JavaScript implementation.

## Features

- Generates L1 Voronoi diagrams using Manhattan distance (|x1-x2| + |y1-y2|)
- Implements divide-and-conquer approach with recursive splitting and merging
- Visualizes both the final diagram and the merge process using matplotlib
- Handles edge cases and square bisectors with point nudging
- Object-oriented design with Site and Bisector classes

## Installation

1. Install the required dependencies:
```bash
pip install -r requirements.txt
```

## Usage

### Quick Demo
Run the demo script to generate a Voronoi diagram with 16 random points:

```bash
python demo.py
```

This will:
1. Generate 16 random points using normal distribution
2. Create the L1 Voronoi diagram
3. Save the visualization to 'manhattan_voronoi.png'
4. Display statistics about the generated cells

### Simple Test
Run a simple test with 4 predefined points:

```bash
python test_simple.py
```

### Custom Usage
```python
from voronoi import generate_l1_voronoi, visualize_voronoi

# Define points
points = [[50, 50], [150, 80], [120, 150], [80, 120]]
width, height = 200, 200

# Generate Voronoi diagram
voronoi_sites = generate_l1_voronoi(points, width, height, nudge_data=True)

# Visualize
visualize_voronoi(voronoi_sites, width, height, merge_process=True)
```

## Algorithm Overview

The implementation follows Lee and Wong's algorithm for L1 Voronoi diagrams:

1. **Point Generation**: Points are generated (or provided) and sorted by x-coordinate
2. **Data Cleaning**: Optional nudging to avoid square bisectors
3. **Recursive Splitting**: The set of points is recursively divided until we have ≤ 2 points
4. **Bisector Generation**: For 2 points, their L1 bisector is computed
5. **Merging**: Sub-diagrams are merged by walking the merge line between them
6. **Polygon Construction**: Voronoi cells are constructed from the bisectors

## Key Functions

- `generate_l1_voronoi()`: Main function to generate the Voronoi diagram
- `recursive_split()`: Handles the divide-and-conquer splitting
- `walk_merge_line()`: Walks the merge line between two sets of sites
- `find_l1_bisector()`: Computes the L1 bisector between two sites
- `visualize_voronoi()`: Creates matplotlib visualizations

## Implementation Notes

- The algorithm uses Manhattan (L1) distance instead of Euclidean distance
- Square bisectors occur when two points are at equal L1 distance in both x and y dimensions
- The `nudge_data` parameter slightly adjusts points to avoid square bisectors
- Merge process visualization shows the recursive construction with color-coded levels
- For stability, start with fewer points (8-16) and increase as needed

## Known Limitations

- May have issues with certain point configurations that create square bisectors
- Some polygon points may be duplicated (can be cleaned up in post-processing)
- Performance degrades with very large numbers of points (> 100)

## Differences from JavaScript Version

- Uses matplotlib for visualization instead of SVG
- Implements object-oriented design with Site and Bisector classes
- Uses Python's built-in data structures and type hints
- Includes additional error handling and debugging features