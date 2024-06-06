# adapted from code originally written by Dr. Maxie Schmidt
from copy import copy
from typing import List
from functools import cmp_to_key
from numpy import array, mean
from math import sin, cos, atan2, pi
from collections import defaultdict
from constraint import Problem

from svgwrite import Drawing
from svgwrite.path import Path

__major_version__ = "0.0"
__release__ = "1"
__version__ = "%s.%s" % (__major_version__, __release__)

NUMCPUS = 8
FPNUM_DIGITS = 6

#
# Python constants denoting three of the initial tile shapes
#
RHOMB_TILE = 1
TRIANGLE_TILE = 2
SQUARE_TILE = 3
STAR_TILE = 4


def vector(coord: List[float]) -> List[float]:
    return array(coord)


def matrix(input: List[List[float]]) -> List[List[float]]:
    return array(input)


def midpoint(A, B):
    return (A + B) / 2.0


def midpoint2(t, A, B):
    return (1 - float(t)) * A + float(t) * B


def RotationMatrix(theta):
    return matrix([[cos(theta), -1 * sin(theta)], [sin(theta), cos(theta)]])


def round_coordinate(coord):
    return round(coord, FPNUM_DIGITS)


def round_vector(v):
    return V(round_coordinate(X(v)), round_coordinate(Y(v)))


def filter_nonzero(gaplst, round_prec=0.0001):
    return list(filter(lambda gap: abs(gap) > round_prec, gaplst))


def add(A: vector, B: vector) -> vector:
    return vector([A[0] + B[0], A[1] + B[1]])


def sub(A: vector, B: vector) -> vector:
    return vector([A[0] - B[0], A[1] - B[1]])


import numpy as np


def unit_vector(vector):
    """Returns the unit vector of the vector."""
    return vector / np.linalg.norm(vector)


def angle(A, B, origin):
    norm_a = sub(A, origin)
    norm_b = sub(B, origin)
    v1_u = unit_vector(norm_a)
    v2_u = unit_vector(norm_b)
    x_axis = unit_vector(vector([1, 0]))
    angle1 = np.arccos(np.clip(np.dot(v1_u, x_axis), -1.0, 1.0))

    if norm_a[1] < 0:
        angle1 *= -1
    angle2 = np.arccos(np.clip(np.dot(v2_u, x_axis), -1.0, 1.0))
    if norm_b[1] < 0:
        angle2 *= -1
    diff_angle = angle2 - angle1

    if diff_angle < -pi:
        diff_angle += 2 * pi
    if diff_angle > pi:
        diff_angle -= 2 * pi
    return diff_angle


def project(
    endpoint: vector, startpoint: vector, angle: float, length: float
) -> vector:
    # project a point based on a vector and an angle and length
    beta = atan2((endpoint[1] - startpoint[1]), endpoint[0] - startpoint[0])
    new_angle = angle + beta
    x1 = length * cos(new_angle)
    y1 = length * sin(new_angle)
    return vector([startpoint[0] + x1, startpoint[1] + y1])


def get_solutionsXY(solns):
    if solns == None:
        return []
    XYsol_func = lambda i: [solns[i][0].rhs().real, solns[i][1].rhs().real]
    XYsols = map(XYsol_func, range(0, len(solns)))
    return XYsols


## edist
# Computes the Euclidean distance between two points
# @param p0    The first 2D point
# @param p1    The second 2D point
# @param sqpow The optional power of the distance function
#              (defaults to the normal sqpow = 1/2)
# @return      The (power modified) Euclidean distance between the points
##
def edist(p0, p1=[0, 0], sqpow=0.5):
    (x1, y1, x2, y2) = (p0[0], p0[1], p1[0], p1[1])
    return (((x1 - x2) ** 2) + ((y1 - y2) ** 2)) ** sqpow


##def


## sort_points_complex
# Sorts a list of 2D points, or tuples, using the numpy sort_complex routine
# @param points_list A list of 2D points
# @return            A sorted list containing the original points
#                    sorted first with respect to the first coordinates, then
#                    with respect to the second
##
def sort_points_complex(points_list):
    complex_points = []
    for p in points_list:
        x, y = p[0], p[1]
        cp = np.complex(x, y)
        complex_points.append(cp)

    complex_points = list(np.sort_complex(complex_points))

    points_list = []
    for cp in complex_points:
        re, im = cp.real, cp.imag
        points_list.append(vector([re, im]))

    return points_list


## sort_points_1D
# Sorts a list of one-dimensional (i.e., no tuples) elements
# @param points_list     A list of 1D elements to be sorted
# @param sort_by_ycoords Always ignored
# @return                A sorted list of the elements
##
def sort_points_1D(points_list, sort_by_ycoords=False):
    return list(np.sort(points_list))


## unique_points
# Computes a list of distinct tuples
# @param points_list  A list of pairs
# @param perform_sort Indicates whether to sort the list of points before the
#                     search for distinct pairs
#                     (note the list must be pre-sorted for this function
#                      to correctly determine the unique points in the list)
# @return             A list of the distinct tuples in the list
##
def unique_points(points_list, perform_sort=True):
    if perform_sort == True:
        # points_list = sort_points_complex(points_list)
        points_list = sorted(points_list)
    ## if

    ## Round floating point numbers to avoid zero gaps from non-distinct
    ## points that are numerically very close:
    points_list = map(round_vector, points_list)

    (last_xc, last_yc) = points_list[0]
    uidx = 1
    while uidx < len(points_list):
        (px, py) = points_list[uidx]
        if last_xc != px or last_yc != py:
            (last_xc, last_yc) = (px, py)
            uidx += 1
        else:
            points_list.pop(uidx)
            ## if
    ## while

    return points_list


## unique_points_1D
# Determines the unique 1D elements contained in the list
# @param points_list  A list of 1D (i.e., no tuples) elements
# @param perform_sort Indicates whether to sort the list of points before the
#                     search for distinct elements
#                     (note the list must be pre-sorted for this function
#                      to correctly determine the unique points in the list)
# @return             A list of the distinct elements in the list
##
def unique_points_1D(points_list, perform_sort=True):
    if perform_sort == True:
        points_list = sort_points_1D(points_list)
        ## if

    ## Round floating point numbers to avoid zero gaps from non-distinct
    ## points that are numerically very close:
    points_list = map(round_coordinate, points_list)

    last_point = points_list[0]
    uidx = 1
    while uidx < len(points_list):
        point = points_list[uidx]
        if last_point != point:
            last_point = point
            uidx += 1
        else:
            points_list.pop(uidx)
            ## if
    ## while

    return points_list


## Tiling
# A super class intended to generate derived individual tiling classes
# implemented in the program
##
class Tiling(object):
    ## __init__
    # The initialization function for the Tiling class
    # @param num_steps_N     Parameter number of substitution or replacement
    #                        steps used in generating the sub-tiling
    # @param tiling_name_str A unique string identifier for the tiling
    ##
    def __init__(self, num_steps_N, tiling_name_str):
        self.num_steps = num_steps_N
        self.tiling_name = tiling_name_str
        self.tiles = []
        self.draw_debug = False

    ## N
    # Provides access to the num_steps property of the tiling
    ##
    @property
    def N(self):
        return copy(self.num_steps)

        ## name

    # Provides access to the tiling_name property of the tiling
    ##
    @property
    def name(self):
        return copy(self.tiling_name)

        ## desc

    # Returns a string description of the tiling
    # (intended to be overridden by the sub-classes)
    ##
    def desc(self):
        return "<Tiling Description>"

        ## get_tiles

    # Returns the only tile in the tiling after zero steps
    #
    def get_initial_tile(self):
        return self.INIT_TILE

        # get_next_tiling

    # Gets the next list of tiles after one subsequent substitution step
    # @param prev_tiles        A list of the tiles after one step back
    # @return                  A list of tiles after one more step
    #
    def get_next_tiling(self, prev_tiles):
        next_tiles = []
        for tile in prev_tiles:
            subtiles = tile.to_subtiles()
            next_tiles.extend(subtiles)

        return next_tiles

    # get_tiles
    # Gets the polygonal Ammann tiles after N steps
    # @return A list of tiles in the computed substitution tiling
    #
    def get_tiles(self):
        tile_list = self.get_initial_tile()
        for n in range(1, self.N):
            next_tiles_list = self.get_next_tiling(tile_list)
            tile_list = next_tiles_list
        # TODO: figure out how to combine the triangle tiles
        rtiles_list = []
        for idx, atile in enumerate(tile_list):
            rtiles_list.append(atile)

        return rtiles_list

    def gen_svg(self):
        if not self.tiles:
            return
        size = 200
        [x_min, x_max, y_min, y_max] = get_bounding_box(self.tiles, size)
        view_box = [-x_min, y_min, x_min * 2, y_max - y_min]
        squares = [tile for tile in self.tiles if tile.tile_type == SQUARE_TILE]
        rhomboids = [tile for tile in self.tiles if tile.tile_type == RHOMB_TILE]
        triangles = [tile for tile in self.tiles if tile.tile_type == TRIANGLE_TILE]
        dwg = Drawing(
            filename=f"{self.name}_{self.num_steps}.svg",
            viewBox="{} {} {} {}".format(*view_box),
        )

        if self.draw_debug:
            for tile in squares + rhomboids:
                plot_shape(dwg, tile, size, "red")

            counts = defaultdict(int)
            for tile in triangles:
                plot_shape(dwg, tile, size, fill_color="purple")
        else:
            adjacency_matrix = defaultdict(list)
            final_shapes = [
                tile for tile in self.tiles  # if tile.tile_type != TRIANGLE_TILE
            ]
            for i, tile in enumerate(final_shapes):
                tile.id = str(i)

            for i, shape1 in enumerate(final_shapes):
                for j, shape2 in enumerate(final_shapes):
                    if i == j:
                        continue
                    if i in adjacency_matrix[j] or j in adjacency_matrix[i]:
                        continue
                    if are_neighbors(shape1.to_points(), shape2.to_points()):
                        adjacency_matrix[i].append(j)
                        adjacency_matrix[j].append(i)
            num_colors = 4
            colors = [i for i in range(num_colors)]

            problem = Problem()
            for i in range(len(final_shapes)):
                problem.addVariable(f"{i}", colors)
            for i in range(len(final_shapes)):
                for neighbor in adjacency_matrix[i]:
                    problem.addConstraint(lambda a, b: a != b, (str(neighbor), str(i)))
            coloring = problem.getSolution()
            color_names = ["red", "green", "blue", "yellow", "purple", "orange"]

            if not coloring:
                print("no coloring solutions were found!")
                coloring = {}
                for i in range(len(final_shapes)):
                    coloring[str(i)] = i % len(color_names)
            for i, shape in enumerate(final_shapes):
                plot_shape(dwg, shape, size, color_names[coloring[str(i)]])
        dwg.save(pretty=True)
        print(f"there were {len(squares)} squares and {len(rhomboids)} rhomboids")


def plot_shape(dwg, shape, size, fill_color):
    id = shape.id
    points = shape.to_points()
    path_string = "M {} {} ".format(points[0][0] * size, points[0][1] * size)
    for point in points[1:] + points[:1]:
        path_string += "L {} {} ".format(point[0] * size, point[1] * size)
    dwg.add(
        Path(
            d=path_string,
            fill=fill_color,
            id=f"{type(shape).__name__}_{id}",
        )
    )


# sort all the points in the tiles
def xy_sort_function(tile1, tile2):
    minx1 = min(point[0] for point in tile1.to_points())
    minx2 = min(point[0] for point in tile2.to_points())
    miny1 = min(point[1] for point in tile1.to_points())
    miny2 = min(point[1] for point in tile2.to_points())
    if minx1 == minx2:
        return miny2 - miny1
    return minx1 - minx2


def distance(point1, point2):
    return ((point1[0] - point2[0]) ** 2.0 + (point1[1] - point2[1]) ** 2.0) ** 0.5


def get_bounding_box(tiles, size):
    tiles = sorted(tiles, key=cmp_to_key(xy_sort_function))

    x_values = []
    y_values = []
    for tile in tiles:
        for point in tile.to_points():
            x_values.append(point[0] * size)
            y_values.append(point[1] * size)

    y_min = min(y_values)
    y_max = max(y_values)
    x_min = min(x_values)
    x_max = max(x_values)
    return [x_min, x_max, y_min, y_max]


epsilon = 0.05


def are_neighbors(points1, points2):
    sides1 = [(points1[i], points1[i - 1]) for i in range(len(points1))]
    sides2 = [(points2[i], points2[i - 1]) for i in range(len(points2))]
    for side1i, side1 in enumerate(sides1):
        for side2i, side2 in enumerate(sides2):
            distance11 = (
                (side1[0][0] - side2[0][0]) ** 2.0 + (side1[0][1] - side2[0][1]) ** 2.0
            ) ** 0.5
            distance12 = (
                (side1[0][0] - side2[1][0]) ** 2.0 + (side1[0][1] - side2[1][1]) ** 2.0
            ) ** 0.5
            distance21 = (
                (side1[1][0] - side2[0][0]) ** 2.0 + (side1[1][1] - side2[0][1]) ** 2.0
            ) ** 0.5
            distance22 = (
                (side1[1][0] - side2[1][0]) ** 2.0 + (side1[1][1] - side2[1][1]) ** 2.0
            ) ** 0.5
            if distance11 < epsilon / 3 and distance22 < epsilon / 3:
                return True
            elif distance12 < epsilon / 3 and distance21 < epsilon / 3:
                return True


def mean_size(tiles):
    return mean(
        [distance(tile[i - 1], tile[i]) for tile in tiles for i in range(len(tile))]
    )
