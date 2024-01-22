# adapted from code originally written by Dr. Maxie Schmidt
from collections import defaultdict
from copy import copy
from math import sin, cos, sqrt, pi

from functools import cmp_to_key

from tiling import Tiling, vector
from svgwrite import Drawing
from svgwrite.path import Path

## d1, d2, sc, s2p1
# Python floats representing several abbreviations for the constant values
# we need in the routines below
##
d1 = float(1 / (1 + sqrt(2)))
d2 = float(sqrt(2))
sc = float(cos(pi / 4))
s2p1 = float(sqrt(2) + 1)

##
# Python constants denoting three of the initial tile shapes
##
RHOMB_TILE = 1
TRIANGLE_TILE = 2
SQUARE_TILE = 3


## AmmannTile
# A class representing an individual tile in the AmmannOctagon tilings
# The RHOMB_TILE, TRIANGLE_TILE, and SQUARE_TILE variants each have
# different substitution rules implemented as below.
##
class AmmannTile(object):
    ## __init__
    # Initialization function for the AmmannTile class
    # @param x, y, z, w Coordinates for the vertices of the
    #                   rhombic, triangular, or square shaped tile
    ##
    def __init__(self, tile_type, x, y, z, w=vector([0, 0])):
        self.tile_type = tile_type
        self.x = x
        self.y = y
        self.z = z
        self.w = w

    ## to_points
    # Returns a list of the tile's 2D vertex points
    ##
    def to_points(self):
        if self.tile_type == SQUARE_TILE or self.tile_type == RHOMB_TILE:
            return [self.x, self.y, self.z, self.w]
        else:
            return [self.x, self.y, self.z]

    ## to_subtiles_square
    # Returns a list of tiles corresponding to one substitution step
    # forward for this specific tile
    ##
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

    ## to_subtiles_rhomb
    # Returns a list of tiles corresponding to one substitution step
    # forward for this specific tile
    ##
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

    ## to_subtiles_triangle
    # Returns a list of tiles corresponding to one substitution step
    # forward for this specific tile
    ##
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

    ## to_subtiles
    # Returns a list of tiles corresponding to one substitution step
    # forward for this tile (which corresponds to one of the tile shapes)
    # @see AmmannTile.to_subtiles_square
    # @see AmmannTile.to_subtiles_rhomb
    # @see AmmannTile.to_subtiles_triangle
    ##
    def to_subtiles(self):
        if self.tile_type == SQUARE_TILE:
            return self.to_subtiles_square()
        elif self.tile_type == RHOMB_TILE:
            return self.to_subtiles_rhomb()
        else:
            return self.to_subtiles_triangle()

    ## get_initial_rhomb
    # Static method that returns an initially rhombic tile
    ##
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

    ## get_initial_triangle
    # Static method that returns an initially triangular tile
    ##
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

    ## get_initial_square
    # Static method that returns an initially square tile
    ##
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

    ## get_initial_8star
    # Static method that returns an initially 8-star-shaped tile
    ##
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

    ## get_initial_octagon
    # Static method that returns an initially octagonal tile
    ##
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


## Ammann_Tiling
# A Tiling subclass implementing variants of the Ammann-Beenker tiling
# See: http://tilings.math.uni-bielefeld.de/substitution_rules/ammann_beenker_rhomb_triangle
# See: http://demonstrations.wolfram.com/AmmannTiles/
##
class Ammann_Tiling(Tiling):
    ## __init__
    # Initialization function for the GoldenTriangle_Tiling class
    # @param num_steps_N     The number of substitution steps in the tiling
    # @param tiling_name_str Tiling name string
    # @param tiling_type     Should be one of: RHOMB_TILE, TRIANGLE_TILE,
    #                        SQUARE_TILE, 8 (for an 8-star tile), or
    #                        88 (for an octagonal tile)
    ##
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

    ## desc
    # Returns a description of the tiling
    ##
    def desc(self):
        return "Ammann-Beenker tilings for several initial tiles"

        ## get_initial_tile

    # Returns the only tile in the tiling after zero steps
    ##
    def get_initial_tile(self):
        return self.INIT_TILE

        ## get_next_tiling

    # Gets the next list of tiles after one subsequent substitution step
    # @param prev_tiles        A list of the tiles after one step back
    # @return                  A list of tiles after one more step
    ##
    def get_next_tiling(self, prev_tiles):
        next_tiles = []
        for tile in prev_tiles:
            subtiles = tile.to_subtiles()
            next_tiles.extend(subtiles)

        return next_tiles

    ## get_tiles
    # Gets the polygonal Ammann tiles after N steps
    # @return A list of tiles in the computed substitution tiling
    ##
    def get_tiles(self):
        tile_list = self.get_initial_tile()
        for n in range(1, self.N):
            next_tiles_list = self.get_next_tiling(tile_list)
            tile_list = next_tiles_list

        rtiles_list = []
        for idx, atile in enumerate(tile_list):
            rtiles_list.append(atile.to_points())

        return rtiles_list
class Triangle:
    id: int
    points: [array[float]]


if __name__ == "__main__":
    size = 200
    generations = 3
    epsilon = 0.1
    tiles = Ammann_Tiling(generations, "square", SQUARE_TILE).get_tiles()

    # sort all the points in the tiles
    def xy_sort_function(tile1, tile2):
        minx1 = min(point[0] for point in tile1)
        minx2 = min(point[0] for point in tile2)
        miny1 = min(point[1] for point in tile1)
        miny2 = min(point[1] for point in tile2)
        if minx1 == minx2:
            return miny2-miny1
        return minx1-minx2

    tiles = sorted(tiles, key=cmp_to_key(xy_sort_function))
    # sort all the
    x_values = []
    y_values = []
    for tile in tiles:
        for point in tile:
            x_values.append(point[0]*size)
            y_values.append(point[1]*size)
    y_min = min(y_values)
    y_max = max(y_values)
    x_min = min(x_values)
    x_max = max(x_values)
    triangles = [tile for tile in tiles if len(tile) == 3]
    squares = []
    edge_triangles = []
    # merge the triangles
    while triangles:
        triangle = triangles.pop()
        other_triangles = copy(triangles)
        other_i = 0
        while other_triangles:
            other_triangle = other_triangles.pop()
            same_points = []
            corner_points = []
            comparison_points = triangle+other_triangle
            comparison_points = sorted(comparison_points, key=cmp_to_key(lambda x, y: x[0]-y[0] if x[0] != y[0] else -x[1]+y[1]))
            # reduce the points that over top of each other
            only_unique = [comparison_points[0]]
            for i in range(1, len(comparison_points)):
                point = comparison_points[i-1]
                other_point = comparison_points[i]
                _dist = ((point[0] - other_point[0]) ** 2.0 + (point[1] - other_point[1]) ** 2.0) ** 0.5
                if _dist < epsilon:
                    continue
                only_unique.append(other_point)

            if len(only_unique) == 4:
                squares.append(only_unique[:2]+[only_unique[3],only_unique[2]])
                del triangles[other_i]
                break
            other_i += 1
        edge_triangles.append(triangle)

    rhomboids = [tile for tile in tiles if len(tile) == 4]
    view_box = [-x_min, y_min, x_min * 2, y_max - y_min]

    dwg = Drawing(filename=f"amman_tiling_{generations}.svg", viewBox="{} {} {} {}".format(*view_box))

    for tile in rhomboids:
        path_string = "M {} {}".format(tile[0][0]*size, tile[0][1]*size)
        for point in tile[1:] + tile[:1]:
            path_string += "L {} {}".format(point[0]*size, point[1]*size)
        dwg.add(Path(d=path_string, fill="blue", stroke="black", stroke_width=0.1))
    for tile in squares:
        path_string = "M {} {}".format(tile[0][0]*size, tile[0][1]*size)
        for point in tile[1:] + tile[:1]:
            path_string += "L {} {}".format(point[0]*size, point[1]*size)
        dwg.add(Path(d=path_string, fill="red", stroke="black", stroke_width=0.1))

    """
    for tile in edge_triangles:
        path_string = "M {} {}".format(tile[0][0]*size, tile[0][1]*size)
        for point in tile[1:] + tile[:1]:
            path_string += "L {} {}".format(point[0]*size, point[1]*size)
        dwg.add(Path(d=path_string, fill="green", stroke="black", stroke_width=0.1))
    """
    counts = defaultdict(int)
    for tile in tiles:
        id = counts[len(tile)]

        path_string = "M {} {}".format(tile[0][0] * size, tile[0][1] * size)
        for point in tile[1:] + tile[:1]:
            path_string += "L {} {}".format(point[0] * size, point[1] * size)
        dwg.add(Path(d=path_string, fill="none", stroke="purple", stroke_width=0.5, id=f"{len(tile)}_{id}"))
        counts[len(tile)] += 1

    dwg.save(pretty=True)
