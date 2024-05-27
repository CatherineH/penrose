from tiling import project, vector, angle
from math import pi, sin
import pytest


@pytest.mark.parametrize(
    "endpoint,expected_endpoint",
    [
        (vector([1, 0]), vector([sin(pi / 4), sin(pi / 4)])),
        (vector([-1, 0]), vector([-sin(pi / 4), -sin(pi / 4)])),
        (vector([0, 1]), vector([-sin(pi / 4), sin(pi / 4)])),
    ],
)
def test_project(endpoint, expected_endpoint):
    endpoint = project(endpoint, vector([0, 0]), pi / 4, 1)
    assert pytest.approx(endpoint[0]) == expected_endpoint[0]
    assert pytest.approx(endpoint[1]) == expected_endpoint[1]


@pytest.mark.parametrize(
    "endpoint,expected_endpoint",
    [
        (vector([-1, 0]), vector([-sin(pi / 4), sin(pi / 4)])),
    ],
)
def test_project_negative(endpoint, expected_endpoint):
    endpoint = project(endpoint, vector([0, 0]), -pi / 4, 1)
    assert pytest.approx(endpoint[0]) == expected_endpoint[0]
    assert pytest.approx(endpoint[1]) == expected_endpoint[1]


# negative angles are clockwise
@pytest.mark.parametrize(
    "vector1,vector2",
    [
        (vector([-1, 0]), vector([0, 1])),
        (vector([0, 1]), vector([1, 0])),
        (vector([1, 0]), vector([0, -1])),
        (vector([0, -1]), vector([-1, 0])),
    ],
)
def test_square_angles_negative(vector1, vector2):
    angle_between = angle(vector1, vector2, vector([0, 0]))
    assert -pi / 2 == pytest.approx(angle_between)


# positive angles are counterclockwise
@pytest.mark.parametrize(
    "vector1,vector2",
    [
        (vector([0, 1]), vector([-1, 0])),
        (vector([-1, 0]), vector([0, -1])),
        (vector([0, -1]), vector([1, 0])),
        (vector([0, 1]), vector([1, 0])),
    ],
)
def test_square_angles(vector1, vector2):
    angle_between = angle(vector1, vector2, vector([0, 0]))
    assert pi / 2 == pytest.approx(angle_between)
