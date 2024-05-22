from tiling import Tiling, vector, length, add
from math import atan, sin, cos, sqrt, pi


SQUARE_TILE = 1
RHOMB_TILE = 2
STAR_TILE = 3

sc = float(cos(pi / 4))
sin_alpha = float(sin(pi / 4))

class MillarsTile:
    def __init__(self, tile_type, x:vector, y:vector, z:vector, w:vector):
        self.tile_type = tile_type
        self.x = x
        self.y = y
        self.z = z
        self.w = w

    def to_points(self):
        if self.tile_type != STAR_TILE:
            return [self.x, self.y, self.z, self.w]
        else:
            side_size = length(self.x, self.y)
            x1 = side_size/2.0
            y1 = (side_size*sc)/(2.0*sin_alpha)
            midpoint1 = add(self.x, vector(x1, y1))
            midpoint2 = add(self.y, vector(y1, -x1))
            midpoint3 = add(self.z, vector(-x1, -y1))
            midpoint4 = add(self.w, vector(-x1, -y1))
            return [self.x, midpoint1, self.y, midpoint2, self.z, midpoint3, self.w, midpoint4]

    def to_subtiles_square(self):
        return [MillarsTile(STAR_TILE, self.x, self.y, self.z, self.w)]
    
    def to_subtiles_rhomb(self):
        return [MillarsTile(SQUARE_TILE, _, self.y, _, self.w)]