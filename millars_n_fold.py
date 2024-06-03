from tiling import (
    Tiling,
    vector,
    edist,
    midpoint,
    project,
    STAR_TILE,
    SQUARE_TILE,
    RHOMB_TILE,
    TRIANGLE_TILE,
    angle,
    epsilon,
)

from functools import cmp_to_key
from statistics import mean
from math import sin, cos, sqrt, pi
import pytest


sc = float(cos(pi / 4))
sin_alpha = float(sin(pi / 4))
RATIO_SIDE = sqrt(
                (2.0 - sqrt(2.0)) / (2.0 + sqrt(2.0))
            )

TRIANGLE_RATIO_SIDE = sqrt(2.0)*sqrt(2.0-sqrt(2.0))/2.0


class MillarsTile:
    def __init__(self, tile_type, x: vector, y: vector, z: vector, w: vector = None):
        self.tile_type = tile_type
        self.x = x
        self.y = y
        self.z = z
        self.w = w
        self.id = ""
        self.other_triangle_id = None

    def square_check(self):
        center_point = midpoint(self.x, self.z)
        angle_between = angle(self.x, self.y, center_point) / 2.0
        angle_between2 = angle(self.y, self.z, center_point) / 2.0
        angle_between3 = angle(self.z, self.w, center_point) / 2.0
        angle_between4 = angle(self.w, self.x, center_point) / 2.0

        assert (
            pytest.approx(angle_between)
            == pytest.approx(angle_between2)
            == pytest.approx(angle_between3)
            == pytest.approx(angle_between4)
        ), f"angles not equal! {angle_between} {angle_between2} {angle_between3} {angle_between4} {self.x} {self.y} {self.z} {self.w}"
        return angle_between

    def rhombus_check(self):
        center_point = midpoint(self.x, self.z)
        angles = [0] * 4
        angles[0] = angle(self.x, self.y, center_point) / 2.0
        angles[1] = angle(self.y, self.z, center_point) / 2.0
        angles[2] = angle(self.z, self.w, center_point) / 2.0
        angles[3] = angle(self.w, self.x, center_point) / 2.0
        print(f"rhombus check: {angles=}")

        def sign(x):
            return x > 0

        angle_signs = set(sign(_angle) for _angle in angles)
        print(f"rhombus check: {angles=} {angle_signs=}")
        print(
            f"rhombus check 2: {angle(self.y, self.x, self.w)} == {angle(self.w, self.z, self.y)}, {angle(self.x, self.y, self.z)} == {angle(self.z, self.w, self.x)}"
        )
        rhombus_rules = True
        rhombus_rules = rhombus_rules and pytest.approx(
            angle(self.y, self.x, self.w)
        ) == pytest.approx(angle(self.w, self.z, self.y))
        rhombus_rules = rhombus_rules and pytest.approx(
            angle(self.x, self.y, self.z)
        ) == pytest.approx(angle(self.z, self.w, self.x))
        return len(angle_signs) == 1 and rhombus_rules

    def to_points(self):
        if self.tile_type == TRIANGLE_TILE:
            return [self.x, self.y, self.z]
        if self.tile_type != STAR_TILE:
            return [self.x, self.y, self.z, self.w]
        else:
            center_point = midpoint(self.x, self.z)
            side = edist(self.x, center_point) * RATIO_SIDE
            # we need to find the sign of the square - i.e. is it clockwise or counterclockwise?
            # if it is clockwise the projection angle must be negative, else positive
            angle_between = self.square_check()
            midpoint1 = project(self.x, center_point, angle_between, side)

            # print("angle between:", angle_between, self.x, side, center_point)
            midpoint2 = project(self.y, center_point, angle_between, side)
            midpoint3 = project(self.z, center_point, angle_between, side)
            midpoint4 = project(self.w, center_point, angle_between, side)
            # print(f"midpoints: {midpoint1} {midpoint2} {midpoint3} {midpoint4}")

            return [
                self.x,
                midpoint1,
                self.y,
                midpoint2,
                self.z,
                midpoint3,
                self.w,
                midpoint4,
            ]

    def average_edge_length(self):
        if self.w is None:
            return mean(
                [edist(self.x, self.y), edist(self.y, self.z), edist(self.z, self.x)]
            )
        else:
            return mean(
                [
                    edist(self.x, self.y),
                    edist(self.y, self.z),
                    edist(self.z, self.w),
                    edist(self.w, self.x),
                ]
            )

    def midpoint(self, p1, inverse=False):
        center_point = midpoint(self.x, self.z)
        side = edist(self.x, center_point) * RATIO_SIDE
        angle_between = self.square_check()
        if inverse:
            angle_between *= -1
        return project(p1, center_point, angle_between, side)

    def to_subtiles_square(self):
        # TODO: triangle edges?
        return [
            MillarsTile(STAR_TILE, self.x, self.y, self.z, self.w),
            MillarsTile(
                TRIANGLE_TILE,
                self.x,
                self.midpoint(self.x),
                self.y,
                self.midpoint(self.x, True),
            ),
            MillarsTile(
                TRIANGLE_TILE,
                self.y,
                self.midpoint(self.y),
                self.z,
                self.midpoint(self.y, True),
            ),
            MillarsTile(
                TRIANGLE_TILE,
                self.z,
                self.midpoint(self.z),
                self.w,
                self.midpoint(self.z, True),
            ),
            MillarsTile(
                TRIANGLE_TILE,
                self.w,
                self.midpoint(self.w),
                self.x,
                self.midpoint(self.w, True),
            ),
        ]

    def to_subtiles_rhomb(self):
        center_point = midpoint(self.y, self.w)
        side = edist(self.y, center_point)
        midpoint1 = project(self.y, center_point, pi / 2, side)
        midpoint2 = project(self.w, center_point, pi / 2, side)

        return [
            MillarsTile(SQUARE_TILE, midpoint2, self.y, midpoint1, self.w),
            MillarsTile(TRIANGLE_TILE, self.x, midpoint2, self.y),
            MillarsTile(TRIANGLE_TILE, self.x, midpoint2, self.z),
            MillarsTile(TRIANGLE_TILE, self.w, midpoint1, self.y),
            MillarsTile(TRIANGLE_TILE, self.w, midpoint1, self.z),
        ]

    def to_subtiles_star(self):
        center_point = midpoint(self.x, self.z)
        points = self.to_points()
        print(f"subtiles star {self.x} {self.y} {self.z} {self.w} {center_point}")
        side_length = edist(center_point, points[1]) * sqrt(2)
        angle_between = angle(self.x, self.y, center_point) / 2

        projected_corner1 = project(points[1], center_point, angle_between, side_length)
        projected_corner2 = project(points[3], center_point, angle_between, side_length)
        projected_corner3 = project(points[5], center_point, angle_between, side_length)
        projected_corner4 = project(points[7], center_point, angle_between, side_length)
        square1 = MillarsTile(
            SQUARE_TILE, center_point, points[3], projected_corner1, points[1]
        )
        square1.square_check()
        square2 = MillarsTile(
            SQUARE_TILE, center_point, points[5], projected_corner2, points[3]
        )
        square2.square_check()
        square3 = MillarsTile(
            SQUARE_TILE, center_point, points[7], projected_corner3, points[5]
        )
        square3.square_check()
        square4 = MillarsTile(
            SQUARE_TILE, center_point, points[1], projected_corner4, points[7]
        )
        square4.square_check()
        return [
            square1,
            square2,
            square3,
            square4,
            MillarsTile(TRIANGLE_TILE, self.x, points[1], projected_corner4),
            MillarsTile(TRIANGLE_TILE, self.x, points[7], projected_corner4),
            MillarsTile(TRIANGLE_TILE, self.y, points[1], projected_corner1),
            MillarsTile(TRIANGLE_TILE, self.y, points[3], projected_corner1),
            MillarsTile(TRIANGLE_TILE, self.z, points[3], projected_corner2),
            MillarsTile(TRIANGLE_TILE, self.z, points[5], projected_corner2),
            MillarsTile(TRIANGLE_TILE, self.w, points[5], projected_corner3),
            MillarsTile(TRIANGLE_TILE, self.w, points[7], projected_corner3),
        ]

    def to_subtiles_triangle(self):
        # three triangles, at thirds along the long side
        # the corner that is not being used is common to all three triangles
        dist_12 = edist(self.x, self.y)
        dist_13 = edist(self.x, self.z)
        dist_23 = edist(self.y, self.z)
        common_point = None
        start_point = None
        end_point = None
        _length = None
        if dist_23 > dist_12 and dist_23 > dist_13:
            common_point = self.x
            start_point = self.y
            end_point = self.z
            _length = dist_12
        elif dist_12 > dist_13 and dist_12 > dist_23:
            common_point = self.z
            start_point = self.x
            end_point = self.y
            _length = dist_13
        else:
            common_point = self.y
            start_point = self.x
            end_point = self.z
            _length = dist_23
        a = project(end_point, start_point, 0, _length * TRIANGLE_RATIO_SIDE)
        b = project(start_point, end_point, 0, _length * TRIANGLE_RATIO_SIDE)
        square_midpoint = midpoint(a, b)
        square_length = edist(a, b)
        square_corner = project(square_midpoint, common_point, 0, square_length)
        return [
            MillarsTile(TRIANGLE_TILE, start_point, common_point, a),
            MillarsTile(SQUARE_TILE, a, common_point, b, square_corner),
            MillarsTile(TRIANGLE_TILE, b, common_point, end_point),
        ]

    def to_subtiles(self):
        if self.tile_type == SQUARE_TILE:
            return self.to_subtiles_square()
        elif self.tile_type == RHOMB_TILE:
            return self.to_subtiles_rhomb()
        elif self.tile_type == STAR_TILE:
            return self.to_subtiles_star()
        elif self.tile_type == TRIANGLE_TILE:
            return self.to_subtiles_triangle()


class MillarsNFoldTiling(Tiling):
    def __init__(self, num_steps_N):
        self.tiling_name = "millars_n_fold"

        super().__init__(num_steps_N, self.tiling_name)
        self.num_steps = num_steps_N
        self.INIT_TILE = [
            MillarsTile(
                STAR_TILE,
                vector([-100, 0]),
                vector([0, 100]),
                vector([100, 0]),
                vector([0, -100]),
            )
        ]

    def get_next_tiling(self, prev_tiles):
        next_tiles = []
        for tile in prev_tiles:
            subtiles = tile.to_subtiles()
            next_tiles.extend(subtiles)
        # merge the triangles into rhomboids
        non_triangles = [tile for tile in next_tiles if tile.tile_type != TRIANGLE_TILE]
        triangles = [tile for tile in next_tiles if tile.tile_type == TRIANGLE_TILE]
        for i, triangle in enumerate(triangles):
            triangle.id = i
        remaining_triangles, rhomboids = merge_triangles(triangles)

        next_tiles = non_triangles + rhomboids + remaining_triangles
        return next_tiles


def merge_triangles(triangles):
    rhomboids = []
    # merge the triangles
    print(f"start {triangles=}")
    auto_increment = 0
    for triangle in triangles:
        if not triangle.id:
            triangle.id = f"generated_id{auto_increment}"
            auto_increment += 1

        if triangle.other_triangle_id:
            print(f"we already know {triangle.other_triangle_id=}")
            continue  # we already know this one
        for other_triangle in triangles:
            if not other_triangle.id:
                other_triangle.id = f"generated_id{auto_increment}"
                auto_increment += 1
            if other_triangle.id == triangle.id:
                continue
            if other_triangle.other_triangle_id:
                continue
            comparison_points = triangle.to_points() + other_triangle.to_points()
            average_size_length = 0.5 * (
                triangle.average_edge_length() + other_triangle.average_edge_length()
            )
            debug = False
            only_unique = triangles_to_square(
                comparison_points, expected_size=average_size_length, debug=debug
            )
            # TODO: we need to determine whether the values are concave, that means they should all be going in the same direction?
            # print(f"only_unique {only_unique}")
            if only_unique:
                x_cand = None
                for only_u in only_unique:
                    if x_cand is None:
                        x_cand = only_u
                    if x_cand[0] > only_u[0]:
                        x_cand = only_u
                only_unique = [
                    entry
                    for entry in only_unique
                    if x_cand[0] != entry[0] or x_cand[1] != entry[1]
                ]
                assert len(only_unique) == 3, f"{only_unique} {x_cand}"
                z_cand = None
                for only_u in only_unique:
                    if z_cand is None:
                        z_cand = only_u
                    if z_cand[0] < only_u[0]:
                        z_cand = only_u
                only_unique = [
                    entry
                    for entry in only_unique
                    if z_cand[0] != entry[0] or z_cand[1] != entry[1]
                ]
                y_cand = None
                for only_u in only_unique:
                    if y_cand is None:
                        y_cand = only_u
                    if y_cand[1] > only_u[1]:
                        y_cand = only_u
                only_unique = [
                    entry
                    for entry in only_unique
                    if y_cand[0] != entry[0] or y_cand[1] != entry[1]
                ]
                w_cand = only_unique[0]
                # we may need to swap xz with yw 
                if edist(x_cand, z_cand) < edist(y_cand, w_cand):
                    x_cand, z_cand, y_cand, w_cand = y_cand, w_cand, x_cand, z_cand
                rhomboid = MillarsTile(
                    RHOMB_TILE,
                    x_cand,
                    y_cand,
                    z_cand,
                    w_cand,
                )
                if rhomboid.rhombus_check():
                    rhomboids.append(rhomboid)
                    triangle.other_triangle_id = other_triangle.id
                    other_triangle.other_triangle_id = triangle.id
                else:
                    print(
                        f"points: {rhomboid.to_points()} were determined to not be a rhombus {triangles=}"
                    )
    remaining_triangles = [tile for tile in triangles if not tile.other_triangle_id]
    return remaining_triangles, rhomboids


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
        return only_unique
    else:
        if debug:
            print(f"only_unique points: {len(only_unique)}")
        return None


if __name__ == "__main__":
    for generations in range(1, 6):
        _tiling = MillarsNFoldTiling(generations)
        _tiling.tiles = _tiling.get_tiles()
        print("tiles", _tiling.tiles)
        _tiling.gen_svg()
