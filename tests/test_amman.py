from numpy import array

from ammann_octogon import triangles_to_square


def test_triangles_to_square():
    comparison_points = [array([0., 0.41421356]), array([0.17157288, 0.58578644]), array([0., 0.58578644]), array([0., 0.41421356]),
     array([0.17157288, 0.58578644]), array([0.17157288, 0.41421356])]

    output = triangles_to_square(comparison_points)
    expected = [
        array([0., 0.58578644]), array([0., 0.41421356]), array(
            [0.17157288, 0.41421356]), array([0.17157288, 0.58578644])]
    for i, point in enumerate(output):
        print("point", point, expected[i])
        assert point[0] == expected[i][0]
        assert point[1] == expected[i][1]