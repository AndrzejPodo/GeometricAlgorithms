import numpy as np

def matrix_3x3(point, point_a, point_b):
    matrix = [[point_a[0], point_a[1],1], [point_b[0], point_b[1],1], [point[0], point[1],1]]
    return matrix

def numpy_det(array):
    return np.linalg.det(np.array(array))

def round_point(p, r=2):
    return (round(p[0],r),round(p[1],r))

def round_line(l, r=2):
    return [round_point(l[0],r), round_point(l[1],r)]


def convert_to_dir_graph(lines):
    if len(lines) == 0:
        return {}
    
    lines = sorted(lines, key = lambda line: line[0][1])
    lines = sorted(lines, key=lambda line:line[0][1], reverse=True)
    
    verticies_adj = {lines[0][0]:{lines[0][1]}}
    
    for line in lines:
        if line[0] in verticies_adj:
            if not(line[1] in verticies_adj[line[0]]):
                verticies_adj[line[0]].add(line[1])
        else:
            verticies_adj[line[0]]={line[1]}
    
    return verticies_adj

def convert_to_undir_graph(lines):
    if len(lines) == 0:
        return {}
    
    lines = sorted(lines, key = lambda line: line[0][1])
    lines = sorted(lines, key=lambda line:line[0][1], reverse=True)
    
    verticies_adj = {lines[0][0]:{lines[0][1]}}
    
    for line in lines:
        if line[0] in verticies_adj:
            verticies_adj[line[0]].add(line[1])
        else:
            verticies_adj[line[0]]={line[1]}
        
        if line[1] in verticies_adj:
            verticies_adj[line[1]].add(line[0])
        else:
            verticies_adj[line[1]]={line[0]}
    
    return verticies_adj

def check_monotonous(lines, paths=False):  
    lines = [round_line(line, r=4) for line in lines]
    lines = [sorted(line, key = lambda point: point[1], reverse = True) for line in lines]
    
    starting_point = max(lines, key = lambda line:line[0][1])[0]
    verticies_adj = convert_to_undir_graph(lines)
    verticies_adj = {k: list(v) for k, v in verticies_adj.items()}

    left = min(verticies_adj[starting_point])
    right = max(verticies_adj[starting_point])
    
    if paths:
        left_path = {left}
        right_path = set()

    while left in verticies_adj and min(verticies_adj[left], key=lambda v: v[1])[1] <= left[1]:
        left = min(verticies_adj[left], key=lambda v: v[1])
        if left in left_path:
            break 
        left_path.add(left)

    if paths and not (right in left_path):
        right_path = {right}
    while right in verticies_adj and right != left and min(verticies_adj[right], key=lambda v: v[1])[1] <= right[1]:
        right = min(verticies_adj[right], key=lambda v: v[1])
        if right in right_path:
            break
        if not (right in left_path):
            right_path.add(right)

    if right == left:
        if paths:
            return (True, left_path, right_path)
        return True
    else:
        return False
    
    
def classify_points(lines):
    lines = [round_line(line, r=4) for line in lines]
    lines = [sorted(line, key = lambda point: point[1], reverse = True) for line in lines]
    
    starting_point = max(lines, key = lambda line:line[0][1])[0]
    verticies_adj = convert_to_undir_graph(lines)
    verticies_adj = {k: list(v) for k, v in verticies_adj.items()}
    
    starting_points = [starting_point]
    end_points = []
    connecting_points = []
    dividing_points = []
    valid_points = []

    next = max(verticies_adj[starting_point])
    prev = starting_point
    mid = next
    next = verticies_adj[mid][0] if verticies_adj[mid][0] != prev else verticies_adj[mid][1]

    while mid != starting_point:
        det = numpy_det(matrix_3x3(next, mid, prev))
        if mid[1] < next[1] and mid[1] < prev[1]:
            if det > 0:
                end_points.append(mid)
            else:
                connecting_points.append(mid)
        elif mid[1] > next[1] and mid[1] > prev[1]:
            if det > 0:
                starting_points.append(mid)
            else:
                dividing_points.append(mid)
        else:
            valid_points.append(mid)

        
        prev = mid
        mid = next
        next = verticies_adj[mid][0] if verticies_adj[mid][0] != prev else verticies_adj[mid][1]
    

    return {'starting':starting_points,'connecting':connecting_points,'dividing':dividing_points, 'valid':valid_points, 'end':end_points}
    
def triangulation(lines):
    lines = [round_line(line, r=4) for line in lines]
    lines = [sorted(line, key = lambda point: point[1], reverse = True) for line in lines]
    
    starting_point = max(lines, key = lambda line:line[0][1])[0]
    monotonous, left_path, right_path = check_monotonous(lines, True)

    all_points = sorted(list(left_path) + list(right_path), key=lambda point: point[1], reverse=True)
    stack = []

    stack.append(starting_point)
    stack.append(all_points[0])
    lines.append([stack[-1], stack[-2]])

    for point in all_points[2:]:
        if (point in left_path and stack[-1] in right_path) or (point in right_path and stack[-1] in left_path):
            first = stack[-1]
            while len(stack) > 0:
                lines.append([point,stack.pop()])
            
            stack.append(max(point, first, key=lambda p: p[1]))
            stack.append(min(point, first, key=lambda p: p[1]))
        else:
            if stack[-1] in left_path:
                u1 = stack.pop()
                u2 = stack.pop()
                while len(stack) > 1 and numpy_det(matrix_3x3(stack[-1], point, stack[-2])) > 0:
                    lines.append([point, u2])
                    if len(stack) == 0:
                        break
                    u1 = u2
                    u2 = stack.pop()
                stack.append(u2)
                stack.append(u1)
                stack.append(point)
            else:
                u1 = stack.pop()
                u2 = stack.pop()
                while len(stack) > 1 and numpy_det(matrix_3x3(stack[-1], point, stack[-2])) < 0:
                    lines.append([point, u2])
                    if len(stack) == 0:
                        break
                    u1 = u2
                    u2 = stack.pop()
                stack.append(u2)
                stack.append(u1)
                stack.append(point)
                # det = numpy_det(matrix_3x3(stack[-1], point, stack[-2]))
                # if det > 0:
                #     stack.append(point)
                # else:
                #     lines.append([stack[-2], point])
                #     stack.pop()
                #     stack.append(point)

    return lines 



#lines = [[(0,0),(1,3)],[(0,0),(2,0)],[(2,0),(1,3)]]
lines = [[(0,0),(1,3)],[(0,0),(2,0)],[(2,5),(1,3)],[(2,0),(2,5)]]
#lines = [[(3,0),(1,3)],[(3,0),(4.1,2)],[(4.1,2),(4,5)],[(4,5),(3.5,4)],[(3.5,4),(2,5)],[(2,5),(1,3)]]
#lines = [[(0.4330657897457, 0.9176604831919952), (0.2689528865198936, 0.6581629341723874)], [(0.2689528865198936, 0.6581629341723874), (0.37096901555215167, 0.5267291106429755)], [(0.37096901555215167, 0.5267291106429755), (0.2733883703908613, 0.4188859733880736)], [(0.2733883703908613, 0.4188859733880736), (0.3443561123263452, 0.3043026400547402)], [(0.3443561123263452, 0.3043026400547402), (0.25342869297150644, 0.15264822829003427)], [(0.25342869297150644, 0.15264822829003427), (0.4375012736166678, 0.06839577730964208)], [(0.4375012736166678, 0.06839577730964208), (0.5616948220037645, 0.21330999299591663)], [(0.5616948220037645, 0.21330999299591663), (0.450807725229571, 0.4222560714272892)], [(0.450807725229571, 0.4222560714272892), (0.5350819187779581, 0.5469496988782696)], [(0.5350819187779581, 0.5469496988782696), (0.4042351445844097, 0.7154546008390541)], [(0.4042351445844097, 0.7154546008390541), (0.47298514458440966, 0.8064472478978775)], [(0.47298514458440966, 0.8064472478978775), (0.4330657897457, 0.9176604831919952)]]


print(triangulation(lines))