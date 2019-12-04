import math
import numpy as np
import kdtree 
from scipy.spatial import distance

EPS = 3

class Graph() :
    def __init__(self):
        self.map = dict()
        self.lines_map = dict()

    def add_triangle(self, triangle):
        self.map[triangle] = []
        if triangle.line_a in self.lines_map:
            t = self.lines_map[triangle.line_a]
        else:
            self.lines_map[triangle.line_a] = []
            t = []
        for tri in t:
            self.map[triangle].append(tri)
            self.map[tri].append(triangle)

        self.lines_map[triangle.line_a].append(triangle)

        if triangle.line_b in self.lines_map:
            t = self.lines_map[triangle.line_b]
        else:
            self.lines_map[triangle.line_b] = []
            t = []
        for tri in t:
            self.map[triangle].append(tri)
            self.map[tri].append(triangle)
            
        self.lines_map[triangle.line_b].append(triangle)

        if triangle.line_c in self.lines_map:
            t = self.lines_map[triangle.line_c]
        else:
            self.lines_map[triangle.line_c] = []
            t = []
        for tri in t:
            self.map[triangle].append(tri)
            self.map[tri].append(triangle)
            
        self.lines_map[triangle.line_c].append(triangle)

    def remove_triangle(self, triangle):
        if(triangle in self.map):
            self.lines_map[triangle.line_a].remove(triangle)
            self.lines_map[triangle.line_b].remove(triangle)
            self.lines_map[triangle.line_c].remove(triangle)

            for t in self.get_neighbours(triangle):
                self.map[t].remove(triangle)
            
            del self.map[triangle]

    def get_neighbours(self, triangle):
        if(triangle not in self.map):
             return []
        return self.map[triangle]

    def to_list(self):
        return self.map.keys()

class TriangleIterator:
    def __init__(self, triangle_graph, triangle, point):
        self.triangle_graph = triangle_graph
        self._triangle = triangle
        self.queue = [triangle]
        self.visited = {triangle}
        self.point = point

    def __next__(self):
        if len(self.queue) > 0:
            triangle_to_ret = self.queue.pop(0)
            neighbours = self.triangle_graph.get_neighbours(triangle_to_ret)
    
            bad_neighbours = list(filter(lambda tri: tri not in self.visited and tri.in_circumcil(self.point),neighbours))
            for tri in bad_neighbours:
                self.visited.add(tri)
            self.queue = self.queue + bad_neighbours
            return triangle_to_ret
        else:
            raise StopIteration


class TrianglesSet:
    def __init__(self):
        self.graph = Graph()
        self.centerToTriangleMapping = dict()
        self.kdtree = kdtree.create(dimensions=2)
        self.triangle = None
        self.point = None
        
    def __iter__(self):
        return TriangleIterator(self.graph, self.triangle, self.point)

    def update(self, point):
        removed_points = []
        
        tri_center = self.kdtree.search_nn(point.as_tuple())[0].data
        self.triangle = self.centerToTriangleMapping[tri_center]
        
        while not self.triangle.in_circumcil(point):
            self.kdtree = self.kdtree.remove(tri_center)
            removed_points.append(tri_center)
            debug_var = self.kdtree.search_nn(point.as_tuple())
            tri_center = debug_var[0].data
            #kdtree.visualize(self.kdtree)
            self.triangle = self.centerToTriangleMapping[tri_center]
        
        for p in removed_points:
            self.kdtree.add(p)
            
        self.point = point

    def add(self, triangle):
        self.centerToTriangleMapping[triangle.circle.center.as_tuple()] = triangle
        self.kdtree.add(triangle.circle.center.as_tuple())
        self.graph.add_triangle(triangle)

    def remove(self, triangle):
        del self.centerToTriangleMapping[triangle.circle.center.as_tuple()]
        self.kdtree = self.kdtree.remove(triangle.circle.center.as_tuple())
        self.graph.remove_triangle(triangle)
        # remove also from quad tree
    def to_list(self):
        return self.graph.to_list()

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
    
    def __eq__(self, value):
        return self.as_tuple == value.as_tuple()

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
        return hash(self.as_raw())
    
    def __str__(self):
        return str([self.point_a.as_tuple(), self.point_b.as_tuple()])
    def __eq__(self, value):
        return self.as_raw() == value.as_raw()

    def as_raw(self):
        return (self.point_a.as_tuple(), self.point_b.as_tuple())

class Circle:
    def __init__(self, r, c):
        self.radius = r
        self.center = c
        self.center.round_cords(EPS)

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
    
    def __eq__(self, value):
        return self.circle.center.as_tuple() == value.circle.center.as_tuple() 

def delunay(points):
    #triangulation = set()

    triangulation = TrianglesSet()

    points = [Point(point[0], point[1]) for point in points]
    
    for point in points:
        point.round_cords(EPS)
        
    max_x = max(points, key=lambda point: point.x).x+1
    min_x = min(points, key=lambda point: point.x).x-1
    max_y = max(points, key=lambda point: point.y).y+1
    min_y = min(points, key=lambda point: point.y).y-1

    v_point_1 = Point(round(min_x-(max_x-min_x)/2,EPS), min_y)
    v_point_2 = Point(round(max_x+(max_x-min_x)/2,EPS), min_y)
    v_point_3 = Point(round((max_x+min_x)/2,EPS), round(2*max_y-min_y, EPS))
    # v_point_4 = Point(min_x-1, max_y+1)

    bounding_verticies = {v_point_1.as_tuple(), v_point_2.as_tuple(), v_point_3.as_tuple()}

    points = points + [v_point_1, v_point_2, v_point_3]

    super_tri_1 = Triangle(v_point_1, v_point_2, v_point_3)
    # super_tri_2 = Triangle(v_point_2, v_point_3, v_point_4)

    triangulation.add(super_tri_1)
    # triangulation.add(super_tri_2)

    for point in points:
        polygon = []
        bad_triangles = []
        triangulation.update(point)
        for tri in triangulation:
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

    for tri in triangulation.to_list():
        if tri.point_a.as_tuple() in bounding_verticies or tri.point_b.as_tuple() in bounding_verticies or tri.point_c.as_tuple() in bounding_verticies:
            pass
        else:
            result.append(tri)


    return result


# points = [(94.0763127264591, 28.428745676369658), (-46.497128693305065, -37.40032718707895), (47.400719116640744, 21.819910803485527), (-49.22757698069056, 66.50788135799692), (-16.68004719623839, -67.4448804280598), (-80.35560829929767, -23.800561071483344)]
# points = [(76.2531699979983, -98.05092993300313), (82.24739551803415, 34.69823657711024), (86.69393854638233, -49.28074763752854), (78.26400485407024, 85.43229928985286), (41.19331806407706, 42.23471348335761), (-22.897808718818453, 18.171879521264117)]
#points = [(35.63777803632581, -85.06909870996905), (-48.99820106914661, -1.6016653156986536), (5.881387699744394, 70.71104238180357), (-24.233301365485204, 5.583781768558026), (-59.23286149560904, -19.45005659005838), (64.89026973134989, -34.26736641284134)]
points = [(85.09794612847179, 80.9933710580772), (46.27670000720104, -53.024500903772086), (62.796156161736405, 94.87712415102283), (-7.180721786279179, -68.23786434363363), (44.17976003400457, 59.518867168084455), (84.94845393190349, 21.901669400881033)]
print(len(delunay(points)))