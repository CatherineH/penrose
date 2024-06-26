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
RATIO_SIDE = sqrt((2.0 - sqrt(2.0)) / (2.0 + sqrt(2.0)))

TRIANGLE_RATIO_SIDE = sqrt(2.0) * sqrt(2.0 - sqrt(2.0)) / 2.0


class MillarsTile:
    def __init__(self, tile_type, x: vector, y: vector, z: vector, w: vector = None):
        self.tile_type = tile_type
        self.x = x
        self.y = y
        self.z = z
        self.w = w
        self.id = ""
        self.other_triangle_id = None

    def star_triangle_check(self):
        angles = [0] * 3
        angles[0] = abs(angle(self.x, self.z, self.y))
        angles[1] = abs(angle(self.y, self.x, self.z))
        angles[2] = abs(angle(self.z, self.y, self.x))
        angles.sort()
        assert (
            pytest.approx(angles[0]) == pi / 8.0
        ), f"angles don't look right! {angles[0]}"
        assert (
            pytest.approx(angles[1]) == pi / 8.0
        ), f"angles don't look right! {angles}"
        assert (
            pytest.approx(angles[2]) == 3 * pi / 4.0
        ), f"angles don't look right! {angles[2]}"

    def square_check(self):
        center_point = midpoint(self.x, self.z)
        angle_between = angle(self.x, self.y, center_point) / 2.0
        angle_between2 = angle(self.y, self.z, center_point) / 2.0
        angle_between3 = angle(self.z, self.w, center_point) / 2.0
        angle_between4 = angle(self.w, self.x, center_point) / 2.0
        """
        assert (
            pytest.approx(edist(self.x, self.y))
            == pytest.approx(edist(self.y, self.z))
            == pytest.approx(edist(self.z, self.w))
            == pytest.approx(edist(self.w, self.x))
        )
        """
        if not (
            pytest.approx(angle_between)
            == pytest.approx(angle_between2)
            == pytest.approx(angle_between3)
            == pytest.approx(angle_between4)
        ):
            print(
                f"angles not equal! {angle_between} {angle_between2} {angle_between3} {angle_between4} {self.x} {self.y} {self.z} {self.w}"
            )
            return None
        return angle_between

    def rhombus_check(self):
        center_point = midpoint(self.x, self.z)
        angles = [0] * 4
        angles[0] = angle(self.x, self.y, center_point) / 2.0
        angles[1] = angle(self.y, self.z, center_point) / 2.0
        angles[2] = angle(self.z, self.w, center_point) / 2.0
        angles[3] = angle(self.w, self.x, center_point) / 2.0

        def sign(x):
            return x > 0

        angle_signs = set(sign(_angle) for _angle in angles)
        corner_angles = [
            angle(self.w, self.y, self.x),
            angle(self.w, self.z, self.y),
            angle(self.y, self.w, self.z),
            angle(self.z, self.y, self.w),
        ]
        corner_angles = sorted([abs(corner_angle) for corner_angle in corner_angles])
        rhombus_rules = True
        rhombus_rules = rhombus_rules and pytest.approx(
            corner_angles[0]
        ) == pytest.approx(corner_angles[1])
        rhombus_rules = rhombus_rules and pytest.approx(
            corner_angles[2]
        ) == pytest.approx(corner_angles[3])
        rhombus_rules = rhombus_rules and pytest.approx(corner_angles[0]) == pi / 4.0
        # rhombus_rules = (
        #    rhombus_rules and pytest.approx(corner_angles[-1]) == 3.0 * pi / 4.0
        # )

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

            midpoint2 = project(self.y, center_point, angle_between, side)
            midpoint3 = project(self.z, center_point, angle_between, side)
            midpoint4 = project(self.w, center_point, angle_between, side)

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
        # we will be guaranteed from the triangle merge that x and z are the wide sides of the rhombus, but nothing else
        # about the orientation
        center_point = midpoint(self.y, self.w)
        side = min(
            edist(self.y, center_point), edist(self.x, center_point)
        )  # take into account the possible orientations?
        midpoint1 = project(self.y, center_point, pi / 2, side)
        midpoint2 = project(self.w, center_point, pi / 2, side)
        # we want midpoint2 to be closer to self.x than midpoint1
        if edist(self.x, midpoint1) < edist(self.x, midpoint2):
            midpoint1, midpoint2 = midpoint2, midpoint1

        return [
            MillarsTile(SQUARE_TILE, midpoint2, self.y, midpoint1, self.w),
            MillarsTile(TRIANGLE_TILE, self.x, midpoint2, self.y),
            MillarsTile(TRIANGLE_TILE, self.x, midpoint2, self.w),
            MillarsTile(TRIANGLE_TILE, self.z, midpoint1, self.y),
            MillarsTile(TRIANGLE_TILE, self.z, midpoint1, self.w),
        ]

    def to_subtiles_star(self):
        center_point = midpoint(self.x, self.z)
        points = self.to_points()
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
        _square_tile = MillarsTile(SQUARE_TILE, a, common_point, b, square_corner)
        _triangle_angle_1 = abs(angle(start_point, a, common_point))
        _triangle_angle_2 = abs(angle(end_point, b, common_point))

        if not _square_tile.square_check():
            return []
        assert pytest.approx(_triangle_angle_1) == pi / 8.0
        assert pytest.approx(_triangle_angle_2) == pi / 8.0

        return [
            MillarsTile(TRIANGLE_TILE, start_point, common_point, a),
            _square_tile,
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
    auto_increment = 0
    for triangle in triangles:
        if not triangle.id:
            triangle.id = f"generated_id{auto_increment}"
            auto_increment += 1

        if triangle.other_triangle_id:
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
            if only_unique:
                center_point = sum(only_unique) / 4.0
                dists_to_center = []
                for point in only_unique:
                    dists_to_center.append((point, edist(point, center_point)))
                dists_to_center.sort(key=lambda x: x[1])

                x_cand = dists_to_center[2][0]
                z_cand = dists_to_center[3][0]
                y_cand = dists_to_center[0][0]
                w_cand = dists_to_center[1][0]
                # we may need to swap x/z and y/w for orientation ?

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
    for generations in range(1, 8):
        _tiling = MillarsNFoldTiling(generations)
        _tiling.tiles = _tiling.get_tiles()
        _tiling.gen_svg()
