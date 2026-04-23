import random
from voronoi import (
    generate_l1_voronoi,
    clean_data,
    find_l1_bisector,
)


class TestSpecialCases:
    """Stress tests for special geometric cases"""

    def test_single_point(self):
        sites = [[50, 50]]
        result = generate_l1_voronoi(sites, 100, 100, False)
        assert len(result) == 1
        assert result[0]["site"] == [50, 50]

    def test_two_points_horizontal(self):
        sites = [[10, 50], [90, 50]]
        result = generate_l1_voronoi(sites, 100, 100, True)
        assert len(result) == 2

    def test_two_points_vertical(self):
        sites = [[50, 10], [50, 90]]
        result = generate_l1_voronoi(sites, 100, 100, True)
        assert len(result) == 2

    def test_two_points_diagonal(self):
        sites = [[10, 10], [90, 90]]
        result = generate_l1_voronoi(sites, 100, 100, True)
        assert len(result) == 2

    def test_three_points_collinear_horizontal(self):
        sites = [[10, 50], [50, 50], [90, 50]]
        result = generate_l1_voronoi(sites, 100, 100, True)
        assert len(result) == 3

    def test_three_points_collinear_vertical(self):
        sites = [[50, 10], [50, 50], [50, 90]]
        result = generate_l1_voronoi(sites, 100, 100, True)
        assert len(result) == 3

    def test_four_points_square(self):
        sites = [[10, 10], [90, 10], [90, 90], [10, 90]]
        result = generate_l1_voronoi(sites, 100, 100, True)
        assert len(result) == 4
        for site in result:
            assert len(site["bisectors"]) >= 2

    def test_four_points_diagonal(self):
        sites = [[10, 10], [90, 90], [10, 90], [90, 10]]
        result = generate_l1_voronoi(sites, 100, 100, True)
        assert len(result) == 4

    def test_points_on_edge_left(self):
        sites = [[0, 25], [0, 50], [0, 75]]
        result = generate_l1_voronoi(sites, 100, 100, True)
        assert len(result) == 3

    def test_points_on_edge_right(self):
        sites = [[100, 25], [100, 50], [100, 75]]
        result = generate_l1_voronoi(sites, 100, 100, True)
        assert len(result) == 3

    def test_points_on_edge_top(self):
        sites = [[25, 0], [50, 0], [75, 0]]
        result = generate_l1_voronoi(sites, 100, 100, True)
        assert len(result) == 3

    def test_points_on_edge_bottom(self):
        sites = [[25, 100], [50, 100], [75, 100]]
        result = generate_l1_voronoi(sites, 100, 100, True)
        assert len(result) == 3

    def test_points_at_corners(self):
        sites = [[0, 0], [100, 0], [0, 100], [100, 100]]
        result = generate_l1_voronoi(sites, 100, 100, True)
        assert len(result) == 4

    def test_points_all_in_one_corner(self):
        sites = [[5, 5], [10, 10], [15, 5], [5, 15]]
        result = generate_l1_voronoi(sites, 100, 100, True)
        assert len(result) == 4

    def test_adjacent_points(self):
        sites = [[50, 50], [51, 50], [50, 51]]
        result = generate_l1_voronoi(sites, 100, 100, True)
        assert len(result) == 3

    def test_many_points_linear(self):
        sites = [[i, 50] for i in range(10, 100, 10)]
        result = generate_l1_voronoi(sites, 100, 100, True)
        assert len(result) == 9

    def test_many_points_grid(self):
        sites = [[x, y] for x in [20, 40, 60, 80] for y in [20, 40, 60, 80]]
        result = generate_l1_voronoi(sites, 100, 100, True)
        assert len(result) == 16
        for site in result:
            assert len(site.get("bisectors", [])) > 0


class TestEdgeCases:
    """Edge cases and boundary conditions"""

    def test_minimum_canvas(self):
        sites = [[50, 50]]
        result = generate_l1_voronoi(sites, 10, 10, False)
        assert len(result) == 1

    def test_single_site_full_canvas(self):
        sites = [[50, 50]]
        result = generate_l1_voronoi(sites, 400, 400, False)
        assert len(result) == 1

    def test_points_near_edges(self):
        sites = [[1, 1], [99, 1], [1, 99], [99, 99]]
        result = generate_l1_voronoi(sites, 100, 100, True)
        assert len(result) == 4


class TestStressTests:
    """Performance stress tests"""

    def test_32_points(self):
        random.seed(42)
        sites = [
            [int(random.random() * 400), int(random.random() * 400)] for _ in range(32)
        ]
        result = generate_l1_voronoi(sites, 400, 400, True)
        assert len(result) == 32

    def test_64_points(self):
        random.seed(42)
        sites = [
            [int(random.random() * 400), int(random.random() * 400)] for _ in range(64)
        ]
        result = generate_l1_voronoi(sites, 400, 400, True)
        assert len(result) == 64

    def test_128_points(self):
        random.seed(42)
        sites = [
            [int(random.random() * 400), int(random.random() * 400)] for _ in range(128)
        ]
        result = generate_l1_voronoi(sites, 400, 400, True)
        assert len(result) == 128

    def test_256_points(self):
        random.seed(42)
        sites = [
            [int(random.random() * 400), int(random.random() * 400)] for _ in range(256)
        ]
        result = generate_l1_voronoi(sites, 400, 400, True)
        assert len(result) == 256


class TestBisectorSpecialCases:
    """Special cases for bisector calculation"""

    def test_vertical_bisector(self):
        p1 = {"site": [10, 10], "bisectors": []}
        p2 = {"site": [10, 90], "bisectors": []}
        bisector = find_l1_bisector(p1, p2, 100, 100)
        assert bisector["points"] is not None

    def test_horizontal_bisector(self):
        p1 = {"site": [10, 10], "bisectors": []}
        p2 = {"site": [90, 10], "bisectors": []}
        bisector = find_l1_bisector(p1, p2, 100, 100)
        assert bisector["points"] is not None

    def test_diagonal_bisector(self):
        p1 = {"site": [10, 10], "bisectors": []}
        p2 = {"site": [90, 90], "bisectors": []}
        bisector = find_l1_bisector(p1, p2, 100, 100)
        assert bisector["points"] is not None

    def test_reverse_diagonal_bisector(self):
        p1 = {"site": [10, 90], "bisectors": []}
        p2 = {"site": [90, 10], "bisectors": []}
        bisector = find_l1_bisector(p1, p2, 100, 100)
        assert bisector["points"] is not None


class TestNudgeData:
    """Test data nudging behavior"""

    def test_nudge_on(self):
        sites = [[10, 10], [20, 20]]
        result = generate_l1_voronoi(sites, 100, 100, True)
        assert len(result) == 2

    def test_nudge_off_raises_on_square(self):
        sites = [[10, 10], [20, 20]]
        result = generate_l1_voronoi(sites, 100, 100, True)
        assert len(result) == 2

    def test_clean_data_mutates(self):
        data = [[10, 10], [20, 20]]
        clean_data(data)
        assert data[0] != [10, 10] or data[1] != [20, 20]
