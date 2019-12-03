# import numpy as np
# from scipy.spatial.kdtree import KDTree
# from scipy.spatial import distance

# kdtree = KDTree([(4,4),(0,0),(8,0)])

# x, y = kdtree.query((4.0,3.9))
# d = distance.euclidean(x, y)
# print(d)
# print(kdtree.query_ball_point((4.0,3.9), r=d))

import kdtree 

tree = kdtree.create(dimensions=2)

tree.add((0,0))
tree.add((4,4))
tree.add((8,0))

p = tree.search_nn((7,2))
print(p[0].data)

kdtree.visualize(tree)
tree.remove((4,4))
kdtree.visualize(tree)


