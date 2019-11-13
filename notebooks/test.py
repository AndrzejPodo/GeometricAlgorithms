import math
import numpy as np
import functools

# def cmp_to_key(mycmp):
#     'Convert a cmp= function into a key= function'
#     class K:
#         def __init__(self, obj, *args):
#             self.obj = obj
#         def __lt__(self, other):
#             return mycmp(self.obj, other.obj) < 0
#         def __gt__(self, other):
#             return mycmp(self.obj, other.obj) > 0
#         def __eq__(self, other):
#             return mycmp(self.obj, other.obj) == 0
#         def __le__(self, other):
#             return mycmp(self.obj, other.obj) <= 0
#         def __ge__(self, other):
#             return mycmp(self.obj, other.obj) >= 0
#         def __ne__(self, other):
#             return mycmp(self.obj, other.obj) != 0
#     return K

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


points = [(1,-2),(3,0),(9,2),(2,3),(1,6),(9,6)]

sorted_points = sorted(points, key=lambda point: (point[1],point[0]))
first_point = sorted_points[0]
sorted_points.pop(0)
sorted_points = sorted(sorted_points, key = functools.cmp_to_key(lambda point_a, point_b: classify_point(first_point, point_a, point_b, 0.0001)))
print(first_point)
print(sorted_points)

sorted_points = list(map(lambda point: (point[0], point[1],math.pow(point[1]-first_point[1],2)+math.pow(point[0]-first_point[0],2)), sorted_points))

for i in range(0, len(sorted_points)-1):
    if i < len(sorted_points)-1:
        if classify_point(first_point,(sorted_points[i][0],sorted_points[i][1]),(sorted_points[i+1][0],sorted_points[i+1][1]),0.0001) == 0:
            if(sorted_points[i][2] > sorted_points[i+1][2]):
                del(sorted_points[i+1])
            else:
                del(sorted_points[i])
print(sorted_points)

# print(sorted_points)
# stack = []
# stack.append(first_point)
# stack.append(sorted_points[0])
# stack.append(sorted_points[1])
    
# for point in sorted_points[2:len(sorted_points)]:
#     while classify_point(stack[len(stack)-2], stack[len(stack)-1], point, 0.1) < 0:
#         stack.pop()
#     stack.append(point)

# print(sorted_points)