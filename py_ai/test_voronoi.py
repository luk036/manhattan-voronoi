import pytest
from voronoi import (
    generate_l1_voronoi,
    generate_voronoi_points,
    clean_data,
    distance,
    same_point,
    angle,
    segment_intersection,
)


class TestDistance:
    def test_same_point(self):
        assert distance([0, 0], [0, 0]) == 0

    def test_manhattan_distance(self):
        assert distance([0, 0], [3, 4]) == 7

    def test_manhattan_distance_negative(self):
        assert distance([5, 5], [2, 1]) == 7


class TestSamePoint:
    def test_identical_points(self):
        assert same_point([10, 20], [10, 20]) is True

    def test_different_points(self):
        assert same_point([10, 20], [10, 21]) is False

    def test_different_x(self):
        assert same_point([10, 20], [30, 20]) is False


class TestAngle:
    def test_angle_zero(self):
        a = angle([0, 0], [1, 0])
        assert abs(a) < 0.001

    def test_angle_90_degrees(self):
        a = angle([0, 0], [0, 1])
        assert abs(a - 1.570796) < 0.001

    def test_angle_45_degrees(self):
        a = angle([0, 0], [1, 1])
        assert abs(a - 0.785398) < 0.001


class TestSegmentIntersection:
    def test_intersecting_segments(self):
        l1 = [[0, 0], [10, 10]]
        l2 = [[0, 10], [10, 0]]
        result = segment_intersection(l1, l2)
        assert result is not None

    def test_parallel_segments(self):
        l1 = [[0, 0], [10, 0]]
        l2 = [[0, 5], [10, 5]]
        assert segment_intersection(l1, l2) is None


class TestCleanData:
    def test_clean_data_basic(self):
        data = [[10, 10], [20, 20]]
        result = clean_data(data)
        assert len(result) == 2


class TestGenerateL1Voronoi:
    def test_single_point(self):
        sites = [[10, 10]]
        result = generate_l1_voronoi(sites, 100, 100, False)
        assert len(result) == 1
        assert result[0]["site"] == [10, 10]

    def test_two_points(self):
        sites = [[10, 10], [90, 90]]
        result = generate_l1_voronoi(sites, 100, 100, False)
        assert len(result) == 2
        assert all(s["site"] in [[10, 10], [90, 90]] for s in result)

    def test_four_points(self):
        sites = [[10, 10], [20, 20], [30, 30], [40, 40]]
        result = generate_l1_voronoi(sites, 100, 100, True)
        assert len(result) == 4

    def test_returns_sites_with_bisectors(self):
        sites = [[10, 10], [90, 90]]
        result = generate_l1_voronoi(sites, 100, 100, False)
        assert all("bisectors" in s for s in result)
        assert all(len(s["bisectors"]) > 0 for s in result)

    def test_duplicate_points_raises(self):
        sites = [[10, 10], [10, 10]]
        with pytest.raises((ValueError, ZeroDivisionError)):
            generate_l1_voronoi(sites, 100, 100, False)


class TestGenerateVoronoiPoints:
    def test_basic_callback(self):
        points = [[0, 0], [10, 0]]
        result = generate_voronoi_points(points, 10, 10, distance)
        assert len(result) == 100
