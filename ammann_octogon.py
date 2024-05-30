# adapted from code originally written by Dr. Maxie Schmidt
from collections import defaultdict
from constraint import Problem
from math import atan, sin, cos, sqrt, pi

from functools import cmp_to_key
import sys
from typing import Optional, List, Any


from tiling import (
    Tiling,
    vector,
    get_bounding_box,
    RHOMB_TILE,
    SQUARE_TILE,
    TRIANGLE_TILE,
    mean_size,
    distance,
    epsilon,
)
from svgwrite import Drawing
from svgwrite.path import Path

# d1, d2, sc, s2p1
# Python floats representing several abbreviations for the constant values
# we need in the routines below
#
d1 = float(1 / (1 + sqrt(2)))
d2 = float(sqrt(2))
sc = float(cos(pi / 4))
s2p1 = float(sqrt(2) + 1)


# AmmannTile
# A class representing an individual tile in the AmmannOctagon tilings
# The RHOMB_TILE, TRIANGLE_TILE, and SQUARE_TILE variants each have
# different substitution rules implemented as below.


class AmmannTile(object):
    # Initialization function for the AmmannTile class
    # @param x, y, z, w Coordinates for the vertices of the
    #                   rhombic, triangular, or square shaped tile
    def __init__(self, tile_type, x, y, z, w=vector([0, 0])):
        self.tile_type = tile_type
        self.x = x
        self.y = y
        self.z = z
        self.w = w

    # to_points
    # Returns a list of the tile's 2D vertex points
    def to_points(self):
        if self.tile_type == SQUARE_TILE or self.tile_type == RHOMB_TILE:
            return [self.x, self.y, self.z, self.w]
        else:
            return [self.x, self.y, self.z]

    # Returns a list of tiles corresponding to one substitution step
    # forward for this specific tile
    def to_subtiles_square(self):
        x, y, z, w = self.x, self.y, self.z, self.w
        c1 = x + (d1 * (y - x) + w - x) / d2
        c2 = x + (d1 * ((y - x) + (w - x))) / d2
        c3 = y + (d1 * ((x - y) + (z - y))) / d2
        c4 = d1 * (z - y)
        c5 = d1 * (z - w)
        c6 = d1 * (w - x)
        c7 = d1 * d2 * (y - x)
        return [
            AmmannTile(RHOMB_TILE, x, c2, c1, x + c6),
            AmmannTile(RHOMB_TILE, c1, c1 + c5, z, w + d2 * c5),
            AmmannTile(RHOMB_TILE, y, y + c4, c1 + c5, c3),
            AmmannTile(RHOMB_TILE, c2, x + c7, y, c3),
            AmmannTile(SQUARE_TILE, c1 + c5, c1, c2, c3),
            AmmannTile(TRIANGLE_TILE, x, x + c7, c2),
            AmmannTile(TRIANGLE_TILE, w, x + c6, c1),
            AmmannTile(TRIANGLE_TILE, w, w + d2 * c5, c1),
            AmmannTile(TRIANGLE_TILE, z, y + c4, c1 + c5),
        ]

    # Returns a list of tiles corresponding to one substitution step
    # forward for this specific tile
    def to_subtiles_rhomb(self):
        x, y, z, w = self.x, self.y, self.z, self.w
        c1 = x + d1 * (y - x)
        c2 = d1 * (w - x)
        c3 = d1 * (w - z)
        c4 = d1 * (y - z)
        c5 = d1 * d2 * (z - y)
        c6 = d1 * d2 * (z - w)
        return [
            AmmannTile(RHOMB_TILE, x, c1, c1 + c2, x + c2),
            AmmannTile(RHOMB_TILE, y, z + c3 + c4, w, c1 + c2),
            AmmannTile(RHOMB_TILE, z + c3 + c4, y + c5, z, w + c6),
            AmmannTile(TRIANGLE_TILE, y, c1, c1 + c2),
            AmmannTile(TRIANGLE_TILE, y, y + c5, z + c3 + c4),
            AmmannTile(TRIANGLE_TILE, w, w + c6, z + c3 + c4),
            AmmannTile(TRIANGLE_TILE, w, x + c2, c1 + c2),
        ]

    # Returns a list of tiles corresponding to one substitution step
    # forward for this specific tile
    def to_subtiles_triangle(self):
        x, y, z = self.x, self.y, self.z
        c1 = x + d1 * (y - x) / d2
        c2 = x + (y - x) / d2
        c3 = d1 * (z - x)
        c4 = y + d1 * d2 * (z - y)
        return [
            AmmannTile(RHOMB_TILE, x, c1, c1 + c3, x + c3),
            AmmannTile(RHOMB_TILE, c2, c4, z, c1 + c3),
            AmmannTile(TRIANGLE_TILE, c2, c1, c1 + c3),
            AmmannTile(TRIANGLE_TILE, y, c4, c2),
            AmmannTile(TRIANGLE_TILE, z, x + c3, c1 + c3),
        ]

    # Returns a list of tiles corresponding to one substitution step
    # forward for this tile (which corresponds to one of the tile shapes)
    # @see AmmannTile.to_subtiles_square
    # @see AmmannTile.to_subtiles_rhomb
    # @see AmmannTile.to_subtiles_triangle
    def to_subtiles(self):
        if self.tile_type == SQUARE_TILE:
            return self.to_subtiles_square()
        elif self.tile_type == RHOMB_TILE:
            return self.to_subtiles_rhomb()
        else:
            return self.to_subtiles_triangle()

    # Static method that returns an initially rhombic tile
    @staticmethod
    def get_initial_rhomb():
        return [
            AmmannTile(
                RHOMB_TILE,
                vector([0, 0]),
                vector([1, 0]),
                vector([1 + float(cos(pi / 4)), float(sin(pi / 4))]),
                vector([float(cos(pi / 4)), float(sin(pi / 4))]),
            )
        ]

    # get_initial_triangle
    # Static method that returns an initially triangular tile
    #
    @staticmethod
    def get_initial_triangle():
        return [
            AmmannTile(
                TRIANGLE_TILE,
                vector([0, 0]),
                vector([d2, 0]),
                vector([float(cos(pi / 4)), float(sin(pi / 4))]),
            )
        ]

    # get_initial_square
    # Static method that returns an initially square tile
    #
    @staticmethod
    def get_initial_square():
        return [
            AmmannTile(
                SQUARE_TILE,
                vector([0, 0]),
                vector([1, 0]),
                vector([1, 1]),
                vector([0, 1]),
            )
        ]

    # get_initial_8star
    # Static method that returns an initially 8-star-shaped tile
    #
    @staticmethod
    def get_initial_8star():
        return [
            AmmannTile(
                RHOMB_TILE,
                vector([0, 0]),
                vector([sc, -sc]),
                vector([1 + sc, -sc]),
                vector([1, 0]),
            ),
            AmmannTile(
                RHOMB_TILE,
                vector([0, 0]),
                vector([0, -1]),
                vector([sc, -1 - sc]),
                vector([sc, -sc]),
            ),
            AmmannTile(
                RHOMB_TILE,
                vector([0, 0]),
                vector([-sc, -sc]),
                vector([-sc, -1 - sc]),
                vector([0, -1]),
            ),
            AmmannTile(
                RHOMB_TILE,
                vector([0, 0]),
                vector([-1, 0]),
                vector([-1 - sc, -sc]),
                vector([-sc, -sc]),
            ),
            AmmannTile(
                RHOMB_TILE,
                vector([0, 0]),
                vector([-sc, sc]),
                vector([-1 - sc, sc]),
                vector([-1, 0]),
            ),
            AmmannTile(
                RHOMB_TILE,
                vector([0, 0]),
                vector([0, 1]),
                vector([-sc, 1 + sc]),
                vector([-sc, sc]),
            ),
            AmmannTile(
                RHOMB_TILE,
                vector([0, 0]),
                vector([sc, sc]),
                vector([sc, 1 + sc]),
                vector([0, 1]),
            ),
            AmmannTile(
                RHOMB_TILE,
                vector([0, 0]),
                vector([1, 0]),
                vector([1 + sc, sc]),
                vector([sc, sc]),
            ),
        ]

    # get_initial_octagon
    # Static method that returns an initially octagonal tile
    #
    @staticmethod
    def get_initial_octagon():
        return [
            AmmannTile(
                RHOMB_TILE,
                vector([0, 0]),
                vector([sc, -sc]),
                vector([1 + sc, -sc]),
                vector([1, 0]),
            ),
            AmmannTile(
                RHOMB_TILE,
                vector([0, 0]),
                vector([0, -1]),
                vector([sc, -1 - sc]),
                vector([sc, -sc]),
            ),
            AmmannTile(
                RHOMB_TILE,
                vector([0, 0]),
                vector([-sc, -sc]),
                vector([-sc, -1 - sc]),
                vector([0, -1]),
            ),
            AmmannTile(
                RHOMB_TILE,
                vector([0, 0]),
                vector([-1, 0]),
                vector([-1 - sc, -sc]),
                vector([-sc, -sc]),
            ),
            AmmannTile(
                RHOMB_TILE,
                vector([0, 0]),
                vector([-sc, sc]),
                vector([-1 - sc, sc]),
                vector([-1, 0]),
            ),
            AmmannTile(
                RHOMB_TILE,
                vector([0, 0]),
                vector([0, 1]),
                vector([-sc, 1 + sc]),
                vector([-sc, sc]),
            ),
            AmmannTile(
                RHOMB_TILE,
                vector([0, 0]),
                vector([sc, sc]),
                vector([sc, 1 + sc]),
                vector([0, 1]),
            ),
            AmmannTile(
                RHOMB_TILE,
                vector([0, 0]),
                vector([1, 0]),
                vector([1 + sc, sc]),
                vector([sc, sc]),
            ),
            AmmannTile(
                SQUARE_TILE,
                vector([sc, -1 - sc]),
                vector([1 + sc, -1 - sc]),
                vector([1 + sc, -sc]),
                vector([sc, -sc]),
            ),
            AmmannTile(
                SQUARE_TILE,
                vector([-sc, -1 - sc]),
                vector([0, -s2p1]),
                vector([sc, -1 - sc]),
                vector([0, -1]),
            ),
            AmmannTile(
                SQUARE_TILE,
                vector([-1 - sc, -sc]),
                vector([-1 - sc, -1 - sc]),
                vector([-sc, -1 - sc]),
                vector([-sc, -sc]),
            ),
            AmmannTile(
                SQUARE_TILE,
                vector([-1 - sc, sc]),
                vector([-s2p1, 0]),
                vector([-1 - sc, -sc]),
                vector([-1, 0]),
            ),
            AmmannTile(
                SQUARE_TILE,
                vector([-sc, 1 + sc]),
                vector([-1 - sc, 1 + sc]),
                vector([-1 - sc, sc]),
                vector([-sc, sc]),
            ),
            AmmannTile(
                SQUARE_TILE,
                vector([sc, 1 + sc]),
                vector([0, s2p1]),
                vector([-sc, 1 + sc]),
                vector([0, 1]),
            ),
            AmmannTile(
                SQUARE_TILE,
                vector([1 + sc, sc]),
                vector([1 + sc, 1 + sc]),
                vector([sc, 1 + sc]),
                vector([sc, sc]),
            ),
            AmmannTile(
                SQUARE_TILE,
                vector([1 + sc, -sc]),
                vector([s2p1, 0]),
                vector([1 + sc, sc]),
                vector([1, 0]),
            ),
            AmmannTile(
                RHOMB_TILE,
                vector([1 + sc, -1 - sc]),
                vector([s2p1, -1]),
                vector([s2p1, 0]),
                vector([1 + sc, -sc]),
            ),
            AmmannTile(
                RHOMB_TILE,
                vector([0, -s2p1]),
                vector([1, -s2p1]),
                vector([1 + sc, -1 - sc]),
                vector([sc, -1 - sc]),
            ),
            AmmannTile(
                RHOMB_TILE,
                vector([-1 - sc, -1 - sc]),
                vector([-1, -s2p1]),
                vector([0, -s2p1]),
                vector([-sc, -1 - sc]),
            ),
            AmmannTile(
                RHOMB_TILE,
                vector([-s2p1, 0]),
                vector([-s2p1, -1]),
                vector([-1 - sc, -1 - sc]),
                vector([-1 - sc, -sc]),
            ),
            AmmannTile(
                RHOMB_TILE,
                vector([-1 - sc, 1 + sc]),
                vector([-s2p1, 1]),
                vector([-s2p1, 0]),
                vector([-1 - sc, sc]),
            ),
            AmmannTile(
                RHOMB_TILE,
                vector([0, s2p1]),
                vector([-1, s2p1]),
                vector([-1 - sc, 1 + sc]),
                vector([-sc, 1 + sc]),
            ),
            AmmannTile(
                RHOMB_TILE,
                vector([1 + sc, 1 + sc]),
                vector([1, s2p1]),
                vector([0, s2p1]),
                vector([sc, 1 + sc]),
            ),
            AmmannTile(
                RHOMB_TILE,
                vector([s2p1, 0]),
                vector([s2p1, 1]),
                vector([1 + sc, 1 + sc]),
                vector([1 + sc, sc]),
            ),
        ]


# Ammann_Tiling
# A Tiling subclass implementing variants of the Ammann-Beenker tiling
# See: http://tilings.math.uni-bielefeld.de/substitution_rules/ammann_beenker_rhomb_triangle
# See: http://demonstrations.wolfram.com/AmmannTiles/
#
class Ammann_Tiling(Tiling):
    # __init__
    # Initialization function for the GoldenTriangle_Tiling class
    # @param num_steps_N     The number of substitution steps in the tiling
    # @param tiling_name_str Tiling name string
    # @param tiling_type     Should be one of: RHOMB_TILE, TRIANGLE_TILE,
    #                        SQUARE_TILE, 8 (for an 8-star tile), or
    #                        88 (for an octagonal tile)
    #
    def __init__(self, num_steps_N, tiling_name_str, tiling_type):
        self.num_steps = num_steps_N
        self.tiling_name = tiling_name_str
        self.tiling_type = tiling_type
        if tiling_type == RHOMB_TILE:
            self.INIT_TILE = AmmannTile.get_initial_rhomb()
        elif tiling_type == TRIANGLE_TILE:
            self.INIT_TILE = AmmannTile.get_initial_triangle()
        elif tiling_type == SQUARE_TILE:
            self.INIT_TILE = AmmannTile.get_initial_square()
        elif tiling_type == 8:
            self.INIT_TILE = AmmannTile.get_initial_8star()
        elif tiling_type == 88:
            self.INIT_TILE = AmmannTile.get_initial_octagon()
        else:
            self.INIT_TILE = []

    # desc
    # Returns a description of the tiling
    #
    def desc(self):
        return "Ammann-Beenker tilings for several initial tiles"

        # get_initial_tile


class Triangle:
    id: int
    points: List[Any]
    other_triangle_id: Optional[int]

    def __init__(
        self, id: int, points: List[Any], other_triangle_id: Optional[int] = None
    ):
        self.id = id
        self.points = points
        self.other_triangle_id = other_triangle_id


class Square:
    id: int
    points: List[Any]

    def __init__(self, id: int, points: List[Any]):
        self.id = id
        self.points = points


class Rhombus:
    id: int
    points: List[Any]

    def __init__(self, id: int, points: List[Any]):
        self.id = id
        self.points = points


def is_square(coordinates, debug=False):
    if debug:
        print(f"is_square: {coordinates}")
    distances = []
    for i in range(len(coordinates)):
        point1 = coordinates[i]
        point2 = coordinates[(i + 1) % len(coordinates)]
        distances.append(distance(point1, point2))
    if debug:
        print(f"distances {distances}")
    dx1 = coordinates[0][0] - coordinates[1][0]
    dy1 = coordinates[0][1] - coordinates[1][1]
    dx2 = coordinates[2][0] - coordinates[1][0]
    dy2 = coordinates[2][1] - coordinates[1][1]
    angle = 0
    if dx1 != 0:
        angle += atan(dy1 / dx1)
    else:
        angle += pi / 2.0 if dy1 > 0 else -pi / 2.0
    if dx2 != 0:
        angle -= atan(dy2 / dx2)
    else:
        angle += pi / 2.0 if dy2 > 0 else -pi / 2.0
    if abs(abs(angle) - (pi / 2.0)) > epsilon:
        if debug:
            print(f"angle not 90 degrees {angle}")
        return False
    differences = [
        abs(distances[i] - distances[i - 1])
        for i in range(1, len(distances))
        if abs(distances[i] - distances[i - 1]) > epsilon / 3
    ]
    if debug:
        print(f"is square? {len(differences)} ")
    return not len(differences) > 0


if __name__ == "__main__":
    size = 200
    generations = 5

    tiles = Ammann_Tiling(generations, "square", SQUARE_TILE).get_tiles()

    average_size_length = mean_size(tiles)
    [x_min, x_max, y_min, y_max] = get_bounding_box(tiles, size)

    triangles = [
        Triangle(id=i, points=tile) for i, tile in enumerate(tiles) if len(tile) == 3
    ]
    squares = []
    debug_target = [52, 53]

    rhomboids = [tile for tile in tiles if len(tile) == 4]
    view_box = [-x_min, y_min, x_min * 2, y_max - y_min]
    rhomboids = [Rhombus(id=i, points=rhomboids[i]) for i in range(len(rhomboids))]
    draw_debug = False
    if len(squares) / len(rhomboids) < 0.48:
        print("something went wrong. Generating debug svg")
        draw_debug = True
    dwg = Drawing(
        filename=f"amman_tiling_{generations}.svg",
        viewBox="{} {} {} {}".format(*view_box),
    )

    if draw_debug:
        for tile in squares + rhomboids:
            id = tile.id
            points = tile.points
            path_string = "M {} {} ".format(points[0][0] * size, points[0][1] * size)
            for point in points[1:] + [points[0]]:
                path_string += "L {} {} ".format(point[0] * size, point[1] * size)
            dwg.add(
                Path(
                    d=path_string,
                    fill="red",
                    stroke="black",
                    stroke_width=0.1,
                    id=f"{type(tile).__name__}_{id}",
                )
            )

        counts = defaultdict(int)
        for tile in triangles:
            id = tile.id
            points = tile.points
            path_string = "M {} {} ".format(points[0][0] * size, points[0][1] * size)
            for point in points[1:] + points[:1]:
                path_string += "L {} {} ".format(point[0] * size, point[1] * size)
            dwg.add(
                Path(
                    d=path_string,
                    fill="none",
                    stroke="purple",
                    stroke_width=0.5,
                    id=f"triangle_{id}",
                )
            )
    else:
        adjacency_matrix = defaultdict(list)
        final_shapes = squares + rhomboids
        for i, shape1 in enumerate(final_shapes):
            for j, shape2 in enumerate(final_shapes):
                if i == j:
                    continue
                if i in adjacency_matrix[j] or j in adjacency_matrix[i]:
                    continue
                if are_neighbors(shape1.points, shape2.points):
                    adjacency_matrix[i].append(j)
                    adjacency_matrix[j].append(i)
        num_colors = 4
        colors = [i for i in range(num_colors)]
        print(adjacency_matrix)

        problem = Problem()
        for i in range(len(final_shapes)):
            problem.addVariable(f"{i}", colors)
        for i in range(len(final_shapes)):
            for neighbor in adjacency_matrix[i]:
                problem.addConstraint(lambda a, b: a != b, (str(neighbor), str(i)))
        coloring = problem.getSolution()
        if not coloring:
            print("no coloring solutions were found!")
        color_names = ["red", "green", "blue", "yellow", "purple", "orange"]
        for i, shape in enumerate(final_shapes):
            id = shape.id
            points = shape.points
            path_string = "M {} {} ".format(points[0][0] * size, points[0][1] * size)
            for point in points[1:] + points[:1]:
                path_string += "L {} {} ".format(point[0] * size, point[1] * size)
            dwg.add(
                Path(
                    d=path_string,
                    fill=color_names[coloring[str(i)]],
                    id=f"{type(shape).__name__}_{id}",
                )
            )
    dwg.save(pretty=True)
    print(f"there were {len(squares)} squares and {len(rhomboids)} rhomboids")


def triangles_to_square(comparison_points, expected_size=1.0, debug=False):
    comparison_points = sorted(
        comparison_points,
        key=cmp_to_key(
            lambda x, y: x[0] - y[0]
            if abs(x[0] - y[0]) > epsilon * expected_size
            else -x[1] + y[1]
        ),
    )
    # reduce the points that over top of each other
    only_unique = [comparison_points[0]]
    if debug:
        print(f"expected_size {expected_size}")
        print(f"after sorting: {comparison_points}")
    for i in range(1, len(comparison_points)):
        point = comparison_points[i - 1]
        other_point = comparison_points[i]
        _dist = (
            (point[0] - other_point[0]) ** 2.0 + (point[1] - other_point[1]) ** 2.0
        ) ** 0.5
        if _dist < epsilon * expected_size:
            continue
        else:
            if debug:
                print(f"distance: {_dist}")

        only_unique.append(other_point)
    if len(only_unique) == 4:
        only_unique = only_unique[:2] + [only_unique[3], only_unique[2]]

    else:
        if debug:
            print(f"only_unique points: {len(only_unique)}")
        return None
    if is_square(only_unique, debug):
        return only_unique
    if debug:
        print(f"{only_unique} was not square!")


def merge_triangles(triangles):
    squares = []
    # merge the triangles
    for triangle in triangles:
        if triangle.other_triangle_id:
            continue  # we already know this one
        for other_triangle in triangles:
            if other_triangle.id == triangle.id:
                continue
            if other_triangle.other_triangle_id:
                continue
            same_points = []
            corner_points = []
            comparison_points = triangle.points + other_triangle.points
            debug = False
            if triangle.id in debug_target and other_triangle.id in debug_target:
                print(triangle.id, other_triangle.id, comparison_points)
                debug = True
            only_unique = triangles_to_square(
                comparison_points, expected_size=average_size_length, debug=debug
            )

            if triangle.id in debug_target and other_triangle.id in debug_target:
                if not only_unique:
                    sys.exit()

            if only_unique:
                print(f"matched {triangle.id} to {other_triangle.id}")

                squares.append(
                    Square(id=f"{triangle.id}_{other_triangle.id}", points=only_unique)
                )
                triangle.other_triangle_id = other_triangle.id
                other_triangle.other_triangle_id = triangle.id
                break
    return squares
