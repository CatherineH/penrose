from millars_n_fold import MillarsTile, merge_triangles
from tiling import (
    STAR_TILE,
    vector,
    plot_shape,
    get_bounding_box,
    edist,
    TRIANGLE_TILE,
    RHOMB_TILE,
)
import pytest
from math import sqrt, sin, pi
from svgwrite import Drawing

RATIO = 0.29289321881345254  # 0.611*sin(pi/4)#sqrt((2-sqrt(2))/(2+sqrt(2)))#*sin(pi/4)


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
                vector([-RATIO, -RATIO]),
            ],
        ),
        (
            [
                vector([-1, 0]),
                vector([0, -1]),
                vector([1, 0]),
                vector([0, 1]),
            ],
            [
                vector([-RATIO, -RATIO]),
                vector([RATIO, -RATIO]),
                vector([RATIO, RATIO]),
                vector([-RATIO, RATIO]),
            ],
        ),
    ],
)
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
            assert point == pytest.approx(
                expected_midpoints[i][j]
            ), f"point didn't match: {i} {j}"


@pytest.mark.parametrize(
    "coords",
    [
        (vector([-1, 0]), vector([0, 1]), vector([1, 0]), vector([0, -1])),
        (vector([0, 0]), vector([1, 1]), vector([0, 2]), vector([-1, 1])),
    ],
)
def test_subtiles_star(coords):
    shape = MillarsTile(STAR_TILE, coords[0], coords[1], coords[2], coords[3])
    print(shape.to_points())

    tiles = shape.to_subtiles_star()
    # all of the distances between points in each star should be the same
    prev_edist = None
    for tile in tiles:
        points = tile.to_points()
        for i in range(len(tile.to_points())):
            new_edist = edist(points[i], points[i - 1])
            if prev_edist:
                assert (
                    pytest.approx(prev_edist) == new_edist
                ), f"not a square! {tile.to_points()}"
            prev_edist = new_edist


def test_subtiles_triangle():
    shape = MillarsTile(TRIANGLE_TILE, vector([0, 3]), vector([0, 0]), vector([1, 1.5]))
    triangles = shape.to_subtiles_triangle()

    assert pytest.approx(triangles[0].z[1]) == 2.024344
    assert pytest.approx(triangles[1].z[1]) == 0.975655


@pytest.mark.parametrize(
    "coords",
    [
        (
            vector([7.173814858237198e-13, 11715.7287525381]),
            vector([-1715.7287525380987, 15857.864376269048]),
            vector([0, 20000]),
            vector([1715.7287525380993, 15857.864376269048]),
        ),
    ],
)
def test_subtiles_rhombus(coords):
    rhombus = MillarsTile(
        RHOMB_TILE,
        x=coords[0],
        y=coords[1],
        z=coords[2],
        w=coords[3],
    )
    rhombus.rhombus_check()
    subtiles = rhombus.to_subtiles_rhomb()
    assert subtiles[0].square_check()

    for i in range(1, len(subtiles)):
        print(f"evaluating triangle {i}")
        assert pytest.approx(edist(subtiles[i].x, subtiles[i].y)) == edist(
            subtiles[i].y, subtiles[i].z
        )
        subtiles[i].star_triangle_check()


@pytest.mark.parametrize(
    "triangle1, triangle2",
    [
        (
            MillarsTile(
                TRIANGLE_TILE,
                vector([0, 20000]),
                vector([7.173814858237198e-13, 11715.7287525381]),
                vector([1715.7287525380993, 15857.864376269048]),
            ),
            MillarsTile(
                TRIANGLE_TILE,
                vector([0, 20000]),
                vector([7.173814858237198e-13, 11715.7287525381]),
                vector([-1715.7287525380987, 15857.864376269048]),
            ),
        ),
        (
            MillarsTile(
                TRIANGLE_TILE,
                vector([4142.135623730951, 10000.0]),
                vector([7.173814858237198e-13, 11715.7287525381]),
                vector([5857.86437626905, 5857.86437626905]),
            ),
            MillarsTile(
                TRIANGLE_TILE,
                vector([5857.86437626905, 5857.86437626905]),
                vector([1715.7287525380993, 7573.59312880715]),
                vector([7.173814858237198e-13, 11715.7287525381]),
            ),
        ),
        (
            MillarsTile(
                TRIANGLE_TILE,
                vector([10000.0, 4142.135623730951]),
                vector([10000.0, 1715.7287525380993]),
                vector([11715.7287525381, 0.0]),
            ),
            MillarsTile(
                TRIANGLE_TILE,
                vector([11715.7287525381, 0.0]),
                vector([11715.7287525381, 2426.406871192851]),
                vector([10000.0, 4142.135623730951]),
            ),
        ),
    ],
)
def test_merge_triangles(triangle1, triangle2):
    remaining_triangles, rhombs = merge_triangles([triangle1, triangle2])
    assert remaining_triangles == []
    assert rhombs[0].tile_type == RHOMB_TILE


def test_merge_triangles_no_merge():
    triangle1 = MillarsTile(
        TRIANGLE_TILE, vector([-2, 1]), vector([0, 1]), vector([0, 0])
    )
    triangle1.id = 1
    triangle2 = MillarsTile(
        TRIANGLE_TILE, vector([2, 1]), vector([0, 1]), vector([0, 0])
    )
    triangle2.id = 2
    remaining_triangles, rhombs = merge_triangles([triangle1, triangle2])
    assert remaining_triangles == [triangle1, triangle2]
    assert rhombs == []


"""
def test_rhombus_check():
    rhombus = MillarsTile(
        RHOMB_TILE,
        x=vector([-1, 1]),
        y=vector([0, 0]),
        z=vector([1, 0]),
        w=vector([-1, -1]),
    )
    assert not rhombus.rhombus_check()
    rhombus = MillarsTile(
        RHOMB_TILE,
        x=vector([-1, 1]),
        y=vector([0, 1]),
        z=vector([1, 0]),
        w=vector([0, 0]),
    )
    assert rhombus.rhombus_check()
"""
