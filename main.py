import numpy as np
import math  as math

def numpy_det(array):
    return np.linalg.det(np.array(array))

def generate_matrix_3x3(point, point_a, point_b):
    matrix = [[point_a[0], point_a[1],1], [point_b[0], point_b[1],1], [point[0], point[1],1]]
    return matrix

def classify_point(point, point_a, point_b, epsilon, matrix_generator=generate_matrix_3x3, det=numpy_det):
    determinant = det(matrix_generator(point, point_a, point_b))
    if abs(determinant) > epsilon:
        if determinant < 0:
            return -1
        else:
            return 1
    if abs(determinant) < epsilon:
        return 0

def jarvis(points):
    hull = []
    first_point = max(points)
    hull.append(first_point)
    next_point = None
    
    while next_point != first_point:
        next_point = points[0]
        for point in points:
            if next_point != point and classify_point(point, hull[-1],next_point, 0.001) > 0:
                next_point = point
        hull.append(next_point)
        points.remove(next_point)
    hull.pop()
    return hull

testPoints = [(1,-2),(3,0),(4,3),(6,7),(6,2),(9,2),(9,6),(4,5),(2,3)]

print(jarvis(testPoints))