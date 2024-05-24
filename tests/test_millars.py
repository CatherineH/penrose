from millars_n_fold import MillarsTile
from tiling import STAR_TILE, vector, plot_shape, get_bounding_box
import pytest
from math import sqrt, sin, pi
from svgwrite import Drawing

RATIO = 0.29289321881345254 #0.611*sin(pi/4)#sqrt((2-sqrt(2))/(2+sqrt(2)))#*sin(pi/4)


@pytest.mark.parametrize(
    "coords,expected_midpoints",
    [
        (
            [
               vector([-1, 0]),
               vector([0, 1]),
               vector([1, 0]),
               vector([0, -1]),
            ],
            [
                vector([-RATIO, RATIO]),
                vector([RATIO, RATIO]),
                vector([RATIO, -RATIO]),
                vector([-RATIO, -RATIO])
            ],
        ),])
def test_star_midpoints(coords, expected_midpoints):
    shape = MillarsTile(STAR_TILE, coords[0], coords[1], coords[2], coords[3])
    points = shape.to_points()
    actual_midpoints = [points[1], points[3], points[5], points[7]]
    size = 100
    [x_min, x_max, y_min, y_max] = get_bounding_box([shape], size)
    view_box = [-x_min, y_min, x_min * 2, y_max - y_min]
    dwg = Drawing(
            filename=f"test_star.svg",
            viewBox="{} {} {} {}".format(*view_box),
        )
    plot_shape(dwg, shape, size, "blue")
    dwg.save(pretty=True)
    for i, value in enumerate(actual_midpoints):
        for j, point in enumerate(value):
            assert pytest.approx(expected_midpoints[i][j]) == point, f"point didn't match: {i} {j}"
