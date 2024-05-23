import pytest
from numpy import array

from ammann_octogon import triangles_to_square, is_square, mean_size


@pytest.mark.parametrize(
    "comparison_points,expected",
    [
        (
            [
                array([0.0, 0.41421356]),
                array([0.17157288, 0.58578644]),
                array([0.0, 0.58578644]),
                array([0.0, 0.41421356]),
                array([0.17157288, 0.58578644]),
                array([0.17157288, 0.41421356]),
            ],
            [
                array([0.0, 0.58578644]),
                array([0.0, 0.41421356]),
                array([0.17157288, 0.41421356]),
                array([0.17157288, 0.58578644]),
            ],
        ),
        (
            [
                array([0.0, 0.41421356]),
                array([0.17157288, 0.58578644]),
                array([0.0, 0.58578644]),
                array([0.0, 0.41421356]),
                array([0.17157288, 0.58578644]),
                array([0.17157288, 0.41421356]),
            ],
            [
                array([0.0, 0.58578644]),
                array([0.0, 0.41421356]),
                array([0.17157288, 0.41421356]),
                array([0.17157288, 0.58578644]),
            ],
        ),
        (
            [
                array([0.70710678, 0.53553391]),
                array([0.77817459, 0.46446609]),
                array([0.77817459, 0.53553391]),
                array([0.70710678, 0.53553391]),
                array([0.77817459, 0.46446609]),
                array([0.70710678, 0.46446609]),
            ],
            None,
        ),
        (
            [
                array([0.46446609, 0.29289322]),
                array([0.53553391, 0.36396103]),
                array([0.46446609, 0.36396103]),
                array([0.46446609, 0.29289322]),
                array([0.53553391, 0.36396103]),
                array([0.53553391, 0.29289322]),
            ],
            None,
        ),
        (
            [
                array([0.0, 0.07106781]),
                array([0.02943725, 0.10050506]),
                array([0.02943725, 0.07106781]),
                array([0.0, 0.07106781]),
                array([0.02943725, 0.10050506]),
                array([0.0, 0.10050506]),
            ],
            None,
        ),
    ],
)
def test_triangles_to_square(comparison_points, expected):
    output = triangles_to_square(
        comparison_points, expected_size=mean_size([comparison_points]), debug=True
    )
    assert output
    if expected:
        for i, point in enumerate(output):
            print("point", point, expected[i])
            assert point[0] == expected[i][0]
            assert point[1] == expected[i][1]


@pytest.mark.parametrize(
    "points",
    [
        (
            [
                array([0.17157288, 0.24264069]),
                array([0.17157288, 0.17157288]),
                array([0.17157288, 0.0]),
                array([0.24264069, 0.07106781]),
            ]
        ),
        (
            [
                array([0.05025253, 0.12132034]),
                array([0.05025253, 0.05025253]),
                array([0.17157288, 0.0]),
                array([0.24264069, 0.07106781]),
            ]
        ),
        (
            [
                array([0.87867966, 0.46446609]),
                array([0.75735931, 0.07106781]),
                array([0.92893219, 0.41421356]),
                array([0.82842712, 0.0]),
            ]
        ),
        (
            [
                array([0.80761184, 0.32233047]),
                array([0.82842712, 0.34314575]),
                array([0.85786438, 0.3137085]),
                array([0.80761184, 0.29289322]),
            ]
        ),
    ],
)
def test_not_is_square(points):
    assert not is_square(points, True)


def test_is_square():
    points = [
        array([0.46446609, 0.36396103]),
        array([0.46446609, 0.29289322]),
        array([0.53553391, 0.29289322]),
        array([0.53553391, 0.36396103]),
    ]
    assert is_square(points)
