import math
import numpy as np
import functools
from sortedcontainers import SortedSet
from scipy.stats import linregress
import heapq 

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

def intersect(line_1, line_2, epsilon=0.001):
    line_1 = sorted(line_1)
    line_2 = sorted(line_2)
    slope_1, intercept_1, r_value, p_value, std_err = linregress([line_1[0][0],line_1[1][0]],[line_1[0][1],line_1[1][1]])
    slope_2, intercept_2, r_value, p_value, std_err = linregress([line_2[0][0],line_2[1][0]],[line_2[0][1],line_2[1][1]])
    W = numpy_det([[slope_1, -1],[slope_2, -1]])
    Wx = numpy_det([[-1*intercept_1,-1],[-1*intercept_2,-1]])
    Wy = numpy_det([[slope_1, -1*intercept_1],[slope_2, -1*intercept_2]])
    x = 0
    y = 0
    if(W != 0):
        x = Wx/W
        y = Wy/W
    
    if x < line_1[1][0]+epsilon and x > line_1[0][0]-epsilon and x < line_2[1][0]+epsilon and x > line_2[0][0]-epsilon:
        return (True,(x,y))
    else:
        return (False,(x,y))


class Event:
    def __init__(self, x):
        self.x = x
    def __lt__(self, other):
        return isinstance(other, Event) and self.x < other.x

    
class StartOfLine(Event):
    def __init__(self, x, line):
        super().__init__(x)
        self.line = line

class EndOfLine(Event):
    def __init__(self, x, line):
        super().__init__(x)
        self.line = line
        
class LinesIntersection(Event):
    def __init__(self, x, lines):
        super().__init__(x)
        self.lines = lines
        
class Line:
    def __init__(self, line, id):
        self.line = line
        slope, intercept, r_value, p_value, std_err = linregress([line[0][0],line[1][0]],[line[0][1],line[1][1]])
        self.slope = slope
        self.intercept = intercept
        self.id = id
    def get_y(self, x):
        return self.slope*x+self.intercept

class SwipeState:
    def __init__(self, x):
        self.x = x
    
    def get_x(self):
        return self.x

    def get_previous(self):
        return self.previous

    def set_x(self, x):
        if(abs(x-self.x) > 0.001):
            self.previous = self.x
        self.x = x
        
def swipe(lines):
    result = set([])
    lines = [Line(line, i) for i, line in enumerate(lines)]
    start_of_lines = list(map(lambda line: StartOfLine(line.line[0][0], line), lines))
    end_of_lines = list(map(lambda line: EndOfLine(line.line[1][0], line), lines))
    
    event_heap = []
    intersection_heap = []
    heapq.heapify(event_heap)
    heapq.heapify(intersection_heap)

    for start in start_of_lines:
        heapq.heappush(event_heap, (start.x, start))
    
    for end in end_of_lines:
        heapq.heappush(event_heap, (end.x, end))
    
    state = SwipeState(0)

    def comarator(line_a, line_b):
        result = line_a.get_y(state.get_x())-line_b.get_y(state.get_x())
        if(abs(result)< 0.001):
            result = line_a.slope - line_b.slope
            if(abs(result)< 0.001):
                return 0
            elif result < 0:
                return -1
            else:
                return 1 
        elif result < 0:
            return -1
        else:
            return 1

    state_set = SortedSet(key = functools.cmp_to_key(comarator))
    checked_intersections = set([])

    while event_heap or intersection_heap:
        if len(intersection_heap) > 0 and intersection_heap[0][0] <= event_heap[0][0]+0.0001:
            event = heapq.heappop(intersection_heap)[1]
        else:
            event = heapq.heappop(event_heap)[1]
        

        if isinstance(event, StartOfLine):
            state.set_x(event.x)
            state_set.add(event.line)
            i = state_set.index(event.line)
            if i > 0 and not ((state_set[i].id, state_set[i-1].id) in checked_intersections):
                in_range, point = intersect(state_set[i].line, state_set[i-1].line)
                checked_intersections.add((state_set[i].id, state_set[i-1].id))
                checked_intersections.add((state_set[i-1].id, state_set[i].id))
                if in_range:
                    result.add(point)
                    heapq.heappush(intersection_heap, (point[0], LinesIntersection(point[0], [state_set[i], state_set[i-1]])))
            if i < (len(state_set)-1) and not ((state_set[i].id, state_set[i+1].id) in checked_intersections):
                in_range, point = intersect(state_set[i].line, state_set[i+1].line)
                checked_intersections.add((state_set[i].id, state_set[i+1].id))
                checked_intersections.add((state_set[i+1].id, state_set[i].id))
                if in_range:
                    result.add(point)
                    heapq.heappush(intersection_heap, (point[0], LinesIntersection(point[0], [state_set[i], state_set[i+1]])))        
        
        elif isinstance(event, EndOfLine):
            i = state_set.index(event.line)
            if i > 0 and i < (len(state_set)-1) and not ((state_set[i-1].id, state_set[i+1].id) in checked_intersections):
                in_range, point = intersect(state_set[i-1].line, state_set[i+1].line)
                checked_intersections.add((state_set[i-1].id, state_set[i+1].id))
                checked_intersections.add((state_set[i+1].id, state_set[i-1].id))
                if in_range:
                    result.add(point)
                    heapq.heappush(intersection_heap, (point[0], LinesIntersection(point[0], [state_set[i-1], state_set[i+1]])))

            del state_set[state_set.index(event.line)]
        
        else:
            line_a = event.lines[0]
            line_b = event.lines[1]
            state.set_x(event.x)
            state.set_x(state.get_previous())
            state_set.remove(line_a)
            state_set.remove(line_b)
            state.set_x(state.get_previous())
            state_set.add(line_a)
            state_set.add(line_b)
            idx_a = state_set.index(line_a)
            idx_b = state_set.index(line_b)

            if idx_a < (len(state_set)-1) and not ((state_set[idx_a].id, state_set[idx_a+1].id) in checked_intersections):
                in_range, point = intersect(state_set[idx_a].line, state_set[idx_a+1].line)
                checked_intersections.add((state_set[idx_a].id, state_set[idx_a+1].id))
                checked_intersections.add((state_set[idx_a+1].id, state_set[idx_a].id))
                if in_range and point[0] > state.get_x():
                    result.add(point)
                    heapq.heappush(intersection_heap, (point[0], LinesIntersection(point[0], [state_set[idx_a], state_set[idx_a+1]])))

            if idx_b > 0 and not ((state_set[idx_b].id, state_set[idx_b-1].id) in checked_intersections):
                in_range, point = intersect(state_set[idx_b].line, state_set[idx_b-1].line)
                checked_intersections.add((state_set[idx_b].id, state_set[idx_b-1].id))
                checked_intersections.add((state_set[idx_b-1].id, state_set[idx_b].id))
                if in_range and point[0] > state.get_x():
                    result.add(point)
                    heapq.heappush(intersection_heap, (point[0], LinesIntersection(point[0], [state_set[idx_b], state_set[idx_b-1]])))

            if idx_b < (len(state_set)-1) and not ((state_set[idx_b].id, state_set[idx_b+1].id) in checked_intersections):
                in_range, point = intersect(state_set[idx_b].line, state_set[idx_b+1].line)
                checked_intersections.add((state_set[idx_b].id, state_set[idx_b+1].id))
                checked_intersections.add((state_set[idx_b+1].id, state_set[idx_b].id))
                if in_range and point[0] > state.get_x():
                    result.add(point)
                    heapq.heappush(intersection_heap, (point[0], LinesIntersection(point[0], [state_set[idx_b], state_set[idx_b+1]])))

            if idx_a > 0 and not ((state_set[idx_a].id, state_set[idx_a-1].id) in checked_intersections):
                in_range, point = intersect(state_set[idx_a].line, state_set[idx_a-1].line)
                checked_intersections.add((state_set[idx_a].id, state_set[idx_a-1].id))
                checked_intersections.add((state_set[idx_a-1].id, state_set[idx_a].id))
                if in_range and point[0] > state.get_x():
                    result.add(point)
                    heapq.heappush(intersection_heap, (point[0], LinesIntersection(point[0], [state_set[idx_a], state_set[idx_a-1]])))
        
        print(event.x)
        print(list(map(lambda line: line.line, state_set)))

    return result

#lines = [[(0,0),(3,3)],[(0,3),(3,0)],[(1.5,2),(3,2)]]
#lines = [[(0,2),(7,2)],[(0,1),(7,3)],[(3,0),(7,5)]]
#lines = [[(1,-2),(3,7)],[(0,2),(8,4)],[(1,5),(9,-4)],[(3,3),(11,7)],[(4,5),(6,3)],[(-1,6),(12,6)]]
#lines = [[(0,2),(8,4)],[(1,5),(9,-4)],[(3,3),(11,7)]]
#lines = [[(1,1),(5,2)],[(1,3),(5,5)],[(1,2),(5,1)],[(1,5),(5,3)]]
#lines = [[(0,0),(6,4)],[(4,4),(10,0)],[(-1,2),(5,-2)],[(3,-2),(9,2)],[(3,1),(5,1)]]
lines = [[(4,4),(10,0)],[(3,-2),(9,2)],[(-1,2),(5,-2)],[(0,0),(6,4)],[(3,1),(5,1)]]
print(swipe(lines))