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
    det11 = np.linalg.det(M11)
    if det11 == 0:
        return False, (0,0),(0,0),0
    else:
        M12 = np.delete(np.delete(A,obj=0, axis=0), obj=1, axis = 1)
        M13 = np.delete(np.delete(A,obj=0, axis=0), obj=2, axis = 1)
        M14 = np.delete(np.delete(A,obj=0, axis=0), obj=3, axis = 1)

        x0 = (np.linalg.det(M12))/(2*det11)
        y0 = -(np.linalg.det(M13))/(2*det11)
        r = math.sqrt(math.pow(x0, 2)+math.pow(y0, 2)+np.linalg.det(M14)/det11)

        return True, x0, y0, r

def d(p1, p2):
    return distance.euclidean((p1.x, p1.y),(p2.x, p2.y))

class Point:
    def __init__(self, x, y):
        self.x = x
        self.y = y
    
    def round_cords(self, r=4):
        self.x = round(self.x, r)
        self.y = round(self.y, r)

    def as_tuple(self):
        return (self.x, self.y)
    
    def __str__(self):
        return str(self.as_tuple())
    
    def __hash__(self):
        return hash((self.x, self.y))
    
class Line:
    def __init__(self, a=0,b=0, raw_line = None):
        if raw_line is None:
            self.point_a = min([a,b], key=lambda point: point.as_tuple()) 
            self.point_b = max([a,b], key=lambda point: point.as_tuple()) 
        else:
            self.point_a = Point(raw_line[0][0], raw_line[0][1])
            self.point_b = Point(raw_line[1][0], raw_line[1][1])

    def round_cords(self, r=4):
        self.point_a.round_cords(r)
        self.point_b.round_cords(r)

    def __hash__(self):
        return hash(self.point_a)+hash(self.point_b)
    
    def __str__(self):
        return str([self.point_a.as_tuple(), self.point_b.as_tuple()])

    def as_raw(self):
        return (self.point_a.as_tuple(), self.point_b.as_tuple())

class Circle:
    def __init__(self, r, c):
        self.radius = r
        self.center = c

class Triangle:
    def __init__(self, a, b, c):
        valid, x0, y0, r = circle_on_three_points(a,b,c)
        if not valid:
            raise Exception("Points are collinear")
        self.circle = Circle(r, Point(x0, y0))
        self.point_a = a 
        self.point_b = b 
        self.point_c = c
        self.line_a = Line(a,b)
        self.line_b = Line(b,c)
        self.line_c = Line(c,a)
        self.circle = Circle(r, Point(x0, y0))


    def in_circumcil(self, point):
        if(d(point, self.circle.center) >= self.circle.radius):
            return False
        else:
            return True
    
    def as_raw(self):
        return [self.line_a.as_raw(), self.line_b.as_raw(), self.line_c.as_raw()]

    def __hash__(self):
        return hash((self.point_a.as_tuple(), self.point_b.as_tuple(), self.point_c.as_tuple()))
    
def delunay(points):
    triangulation = set()

    points = [Point(point[0], point[1]) for point in points]
    
    for point in points:
        point.round_cords(5)
        
    max_x = max(points, key=lambda point: point.x).x
    min_x = min(points, key=lambda point: point.x).x
    max_y = max(points, key=lambda point: point.y).y
    min_y = min(points, key=lambda point: point.y).y

    v_point_1 = Point(min_x-1, min_y-1)
    v_point_2 = Point(max_x+1, min_y-1)
    v_point_3 = Point(max_x+1, max_y+1)
    v_point_4 = Point(min_x-1, max_y+1)

    bounding_verticies = {v_point_1.as_tuple(), v_point_2.as_tuple(), v_point_3.as_tuple(), v_point_4.as_tuple()}

    points = points + [v_point_1, v_point_2, v_point_3, v_point_4]

    super_tri_1 = Triangle(v_point_1, v_point_2, v_point_4)
    super_tri_2 = Triangle(v_point_2, v_point_3, v_point_4)

    triangulation.add(super_tri_1)
    triangulation.add(super_tri_2)

    for point in points:
        polygon = []
        bad_triangles = []    
        for tri in list(triangulation):
            if tri.in_circumcil(point): 
                bad_triangles.append(tri)
                if tri.line_a.as_raw() in polygon: ##wiecej niz dwa trojkaty nie beda wspoldzielic krawedzi???
                    polygon.remove(tri.line_a.as_raw())
                else:
                    polygon.append(tri.line_a.as_raw())

                if tri.line_b.as_raw() in polygon:
                    polygon.remove(tri.line_b.as_raw())
                else:
                    polygon.append(tri.line_b.as_raw())

                if tri.line_c.as_raw() in polygon:
                    polygon.remove(tri.line_c.as_raw())
                else:
                    polygon.append(tri.line_c.as_raw())
        
        for tri in bad_triangles:
            triangulation.remove(tri)

        polygon = [Line(raw_line=line) for line in polygon]
        for line in list(polygon):
            try:
                t = Triangle(point, line.point_a, line.point_b)
            except:
                t = None
            
            if not(t is None):
                triangulation.add(t)

    result = []
    for tri in triangulation:
        if tri.point_a.as_tuple() in bounding_verticies or tri.point_b.as_tuple() in bounding_verticies or tri.point_c.as_tuple() in bounding_verticies:
            pass
        else:
            result.append(tri)


    return result




### Test punktu w okregu
# t = Triangle(Point(0,0), Point(0,1), Point(1,0))
# print(hash(t))
# print(t.in_circumcil(Point(0,0)))
# print(t.in_circumcil(Point(0.5,0.5)))
# print(t.in_circumcil(Point(1,1)))

### Test triangulacji
#points = [(0,0),(2,0),(1,1),(1,2)]
points = [(0,0),(2,0),(1,1.75),(1.,2.)]
print(len(delunay(points)))