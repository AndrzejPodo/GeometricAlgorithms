import math
import numpy as np
from scipy.spatial import distance


def circle_on_three_points(p1,p2,p3):
    A = np.array([
        0,                                0,   0,  1, 
        math.pow(p1.x,2)+math.pow(p1.y,2),p1.x,p1.y,1, 
        math.pow(p2.x,2)+math.pow(p2.y,2),p2.x,p2.y,1,
        math.pow(p3.x,2)+math.pow(p3.y,2),p3.x,p3.y,1
    ]).reshape((4,4))

    M11 = np.delete(np.delete(A,obj=0, axis=0), obj=0, axis = 1)
    
    if np.linalg.det(M11) == 0:
        return False
    else:
        M12 = np.delete(np.delete(A,obj=0, axis=0), obj=1, axis = 1)
        M13 = np.delete(np.delete(A,obj=0, axis=0), obj=2, axis = 1)
        M14 = np.delete(np.delete(A,obj=0, axis=0), obj=3, axis = 1)

        x0 = (np.linalg.det(M12))/(2*np.linalg.det(M11))
        y0 = -(np.linalg.det(M13))/(2*np.linalg.det(M11))
        r = math.sqrt(math.pow(x0, 2)+math.pow(y0, 2)+np.linalg.det(M14)/np.linalg.det(M11))

        return x0, y0, r

def d(p1, p2):
    return distance.euclidean((p1.x, p1.y),(p2.x, p2.y))

class Point:
    def __init__(self, x, y):
        self.x = x
        self.y = y
    
    def round_cords(self, r=4):
        self.x = round(self.x, r)
        self.y = round(self.y, r)

class Line:
    def __init__(self, a,b):
        self.point_a = a 
        self.point_b = b

    def round_cords(self, r=4):
        self.point_a.round_cords(r)
        self.point_b.round_cords(r)

class Circle:
    def __init__(self, r, c):
        self.radius = r
        self.center = c

class Triangle:
    def __init__(self, a, b, c):
        self.point_a = a 
        self.point_b = b 
        self.point_c = c
        
        x0, y0, r = circle_on_three_points(a,b,c)
        self.circle = Circle(r, Point(x0, y0))


    def in_circumcil(self, point):
        if(d(point, self.circle.center) >= self.circle.radius):
            return False
        else:
            return True

def delunay(points):
    triangulation = []

    points = [Point(point[0], point[1]) for point in points]

    max_x = max(points, key=lambda point: point.x)[0]
    min_x = min(points, key=lambda point: point.x)[0]
    max_y = max(points, key=lambda point: point.y)[1]
    min_y = min(points, key=lambda point: point.y)[1]

    v_point_1 = Point(min_x, min_y)
    v_point_2 = Point(max_x, min_y)
    v_point_3 = Point(max_x, max_y)
    v_point_4 = Point(min_x, max_y)

    super_tri_1 = Triangle(v_point_1, v_point_2, v_point_4)
    super_tri_2 = Triangle(v_point_2, v_point_3, v_point_4)

    triangulation.append(super_tri_1)
    triangulation.append(super_tri_2)


### Test punktu w okregu
t = Triangle(Point(0,0), Point(0,1), Point(1,0))
print(t.in_circumcil(Point(0,0)))
print(t.in_circumcil(Point(0.5,0.5)))
print(t.in_circumcil(Point(1,1)))