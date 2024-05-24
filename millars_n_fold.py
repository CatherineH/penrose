from tiling import (
    Tiling,
    vector,
    edist,
    add,
    midpoint,
    project,
    STAR_TILE,
    SQUARE_TILE,
    RHOMB_TILE,
)
from math import atan, sin, cos, sqrt, pi, atan2, tan


sc = float(cos(pi / 4))
sin_alpha = float(sin(pi / 4))


class MillarsTile:
    def __init__(self, tile_type, x: vector, y: vector, z: vector, w: vector):
        self.tile_type = tile_type
        self.x = x
        self.y = y
        self.z = z
        self.w = w
        self.id = ""

    def to_points(self):
        if self.tile_type != STAR_TILE:
            return [self.x, self.y, self.z, self.w]
        else:
            center_point = midpoint(self.x, self.z)
            side = edist(self.x, center_point) * sqrt(
                (2.0 - sqrt(2.0)) / (2.0 + sqrt(2.0))
            )
            print(self.x, center_point)

            midpoint1 = project(self.x, center_point, pi / 4, side)
            midpoint2 = project(self.y, center_point, pi / 4, side)
            midpoint3 = project(self.z, center_point, pi / 4, side)
            midpoint4 = project(self.w, center_point, pi / 4, side)

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

    def to_subtiles_square(self):
        # TODO: triangle edges?
        return [MillarsTile(STAR_TILE, self.x, self.y, self.z, self.w)]

    def to_subtiles_rhomb(self):
        center_point = midpoint(self.y, self.w)
        side = edist(self.y, center_point)
        midpoint1 = project(self.y, center_point, pi / 2, side)
        midpoint2 = project(self.w, center_point, pi / 2, side)

        return [MillarsTile(SQUARE_TILE, midpoint2, self.y, midpoint1, self.w)]

    def to_subtiles_star(self):
        center_point = midpoint(self.x, self.z)
        points = self.to_points()
        side_length = edist(center_point, points[1])

        projected_corner1 = project(points[1], center_point, pi / 2, side_length)
        projected_corner2 = project(points[3], center_point, pi / 2, side_length)
        projected_corner3 = project(points[5], center_point, pi / 2, side_length)
        projected_corner4 = project(points[7], center_point, pi / 2, side_length)
        return [
            MillarsTile(
                SQUARE_TILE, center_point, points[3], projected_corner1, points[1]
            ),
            MillarsTile(
                SQUARE_TILE, center_point, points[5], projected_corner2, points[3]
            ),
            MillarsTile(
                SQUARE_TILE, center_point, points[7], projected_corner3, points[5]
            ),
            MillarsTile(
                SQUARE_TILE, center_point, points[1], projected_corner4, points[7]
            ),
        ]

    def to_subtiles(self):
        if self.tile_type == SQUARE_TILE:
            return self.to_subtiles_square()
        elif self.tile_type == RHOMB_TILE:
            return self.to_subtiles_rhomb()
        else:
            return self.to_subtiles_star()


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



if __name__ == "__main__":
    generations = 5
    _tiling = MillarsNFoldTiling(generations)
    _tiling.tiles = _tiling.get_tiles()
    print(_tiling.tiles)
    _tiling.gen_svg()