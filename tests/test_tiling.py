from tiling import project, vector
from math import pi, sin
import pytest

@pytest.mark.parametrize(
    "endpoint,expected_endpoint",
    [
        (vector([1, 0]),vector([sin(pi/4), -sin(pi/4)])),
        (vector([-1, 0]),vector([-sin(pi/4), sin(pi/4)])),
        ])
def test_project(endpoint, expected_endpoint):
    endpoint = project(endpoint, vector([0, 0]), pi/4, 1)
    assert pytest.approx(endpoint[0]) == expected_endpoint[0]
    assert pytest.approx(endpoint[1]) == expected_endpoint[1]