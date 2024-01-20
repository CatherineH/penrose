# make a penrose tiling card
import math

from numpy import average

from shapely.geometry import Polygon

from penrose import PenroseP3, BtileL, psi
from svgwrite import Drawing
from svgwrite.path import Path
from svgpathtools import Line, Path as ptPath, parse_path

# A simple example starting with a BL tile

scale = 100
tiling = PenroseP3(scale, ngen=7)

theta = 2 * math.pi / 5
rot = math.cos(theta) + 1j * math.sin(theta)
A = -scale / 2 + 0j
B = scale / 2 * rot
C = scale / 2 / psi + 0j
tiling.set_initial_tiles([BtileL(A, B, C)])
tiling.make_tiling()
path_string = ""

# this is the margins of the triangles
margins = 0.8

scale_factor = 8

width = 13  # cm
height = 8  # cm

px_to_cm = 35.433


bounds_lines = [
    Line(-(width - height * 1j) * px_to_cm / 2, -(width + height * 1j) * px_to_cm / 2),
    Line(-(width + height * 1j) * px_to_cm / 2, (width - height * 1j) * px_to_cm / 2),
    Line((width - height * 1j) * px_to_cm / 2, (width + height * 1j) * px_to_cm / 2),
    Line((width + height * 1j) * px_to_cm / 2, -(width - height * 1j) * px_to_cm / 2),
]

bounding_box = Polygon([[p.start.real, p.start.imag] for p in bounds_lines])


for element in tiling.elements:
    points = [element.A, element.B, element.C, element.D]
    points = [point * scale_factor for point in points]
    center_point = average(points)
    new_points = [margins * point + (1 - margins) * center_point for point in points]

    element_poly = Polygon([[p.real, p.imag] for p in new_points])
    test_shape = bounding_box.intersection(element_poly)
    if test_shape.is_empty:
        continue

    x, y = test_shape.exterior.coords.xy
    new_points = [x[i] + y[i] * 1j for i in range(len(x))]

    path_string += "M {} {}".format(new_points[0].real, new_points[0].imag)
    for point in new_points[1:] + new_points[:1]:
        path_string += "L {} {}".format(point.real, point.imag)

margins = 1
x_min = width / 2 + margins
y_min = -height * 2 + margins
y_max = y_min + height * 2 + margins * 4

border_format = "M {} {} L {} {} L {} {} L {} {} L {} {}"
border_coords = [
    -x_min,
    y_min,
    x_min,
    y_min,
    x_min,
    y_max,
    -x_min,
    y_max,
    -x_min,
    y_min,
]
y_inner_min = y_min + margins / 2
x_inner_min = x_min - margins / 2
y_inner_max = y_max - margins / 2
inner_border_coords = [
    -x_inner_min,
    y_inner_min,
    x_inner_min,
    y_inner_min,
    x_inner_min,
    y_inner_max,
    -x_inner_min,
    y_inner_max,
    -x_inner_min,
    y_inner_min,
]
border_coords = [p * px_to_cm for p in border_coords]
inner_border_coords = [p * px_to_cm for p in inner_border_coords]

border_string = border_format.format(*border_coords)
inner_border_string = border_format.format(*inner_border_coords)
view_box = [-x_min, y_min, x_min * 2, y_max - y_min]
view_box = [p * px_to_cm for p in view_box]


dwg = Drawing("example1.svg", viewBox="{} {} {} {}".format(*view_box))
cardstock = parse_path(border_string + " " + path_string)
dwg.add(Path(d=inner_border_string, fill="yellow", stroke="none"))
dwg.add(Path(d=ptPath(*cardstock).d(), fill="blue", stroke="none"))
dwg.save(pretty=True)
