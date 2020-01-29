# -*- coding: utf-8 -*-
"""
Code Challenge: Isolated Stars.
(c) CGWA Time Domain Astronomy Group
"""
import heapq
import itertools
import math


class KDTree(object):
    """A Node that contains kd-tree specific data and methods """

    def __init__(self, data=None, left=None, right=None, axis=None,
            sel_axis=None, dimensions=None, index=None):
        """Creates a new node for a kd-tree
        If the node will be used within a tree, the axis and the sel_axis
        function should be supplied.
        sel_axis(axis) is used when creating subnodes of the current node. It
        receives the axis of the parent node and returns the axis of the child
        node."""
        self.data = data
        self.index = index
        self.left = left
        self.right = right
        self.axis = axis
        self.sel_axis = sel_axis
        self.dimensions = dimensions

    def __repr__(self):
        return '<%(cls)s - %(data)s>' % \
            dict(cls=self.__class__.__name__, data=repr(self.data))

    def __nonzero__(self):
        return self.data is not None

    __bool__ = __nonzero__

    def __eq__(self, other):
        if isinstance(other, tuple):
            return self.data == other
        else:
            return self.data == other.data

    def dist(self, point):
        "Euclidean distance between the current Node and the given point"
        r = range(self.dimensions)
        return sum([math.pow(self.data[i] - point[i], 2) for i in r])

    def search_knn(self, point, k):
        """ Return the k nearest neighbors of point and their distances
        point must be an actual point, not a node.
        k is the number of results to return. The actual results can be less
        (if there aren't more nodes to return) or more in case of equal
        distances.
        dist is a distance function, expecting two points and returning a
        distance value. Distance values can be any comparable type.
        The result is an ordered list of (node, distance) tuples.
        """
        results = []
        self._search_node(point, k, results, itertools.count())

        # We sort the final result by the distance in the tuple
        # (<KdTree>, distance).
        return [(node.data, -d) for d, _, node in sorted(results, reverse=True)]


    def _search_node(self, point, k, results, counter):
        if not self:
            return

        #nodeDist = get_dist(self)
        nodeDist = self.dist(point)

        # Add current node to the priority queue if it closer than
        # at least one point in the queue.
        #
        # If the heap is at its capacity, we need to check if the
        # current node is closer than the current farthest node, and if
        # so, replace it.
        item = (-nodeDist, next(counter), self)
        if len(results) >= k:
            if -nodeDist > results[0][0]:
                heapq.heapreplace(results, item)
        else:
            heapq.heappush(results, item)
        # get the splitting plane
        split_plane = self.data[self.axis]
        # get the squared distance between the point and the splitting plane
        # (squared since all distances are squared).
        plane_dist = point[self.axis] - split_plane
        plane_dist2 = plane_dist * plane_dist

        # Search the side of the splitting plane that the point is in
        if point[self.axis] < split_plane:
            if self.left is not None:
                self.left._search_node(point, k, results, counter)
        else:
            if self.right is not None:
                self.right._search_node(point, k, results, counter)

        # Search the other side of the splitting plane if it may contain
        # points closer than the farthest point in the current results.
        if -plane_dist2 > results[0][0] or len(results) < k:
            if point[self.axis] < self.data[self.axis]:
                if self.right is not None:
                    self.right._search_node(point, k, results, counter)
            else:
                if self.left is not None:
                    self.left._search_node(point, k, results, counter)


def create(point_list=None, dimensions=None, axis=0):
    """ Creates a kd-tree from a list of points
    All points in the list must be of the same dimensionality.
    If no point_list is given, an empty tree is created. The number of
    dimensions has to be given instead.
    If both a point_list and dimensions are given, the numbers must agree.
    Axis is the axis on which the root-node should split.
    sel_axis(axis) is used when creating subnodes of a node. It receives the
    axis of the parent node and returns the axis of the child node. """

    dimensions = dimensions or len(point_list[0])

    # by default cycle through the axis
    sel_axis = lambda prev_axis: (prev_axis + 1) % dimensions

    if not point_list:
        return KDTree(sel_axis=sel_axis, axis=axis, dimensions=dimensions)

    # Sort point list and choose median as pivot element
    point_list = list(point_list)
    point_list.sort(key=lambda point: point[axis])
    median = len(point_list) // 2

    loc   = point_list[median]
    left  = create(point_list[:median], dimensions, sel_axis(axis))
    right = create(point_list[median + 1:], dimensions, sel_axis(axis))
    return KDTree(loc, left, right, axis=axis, dimensions=dimensions)

def isolated(star_list, dist_min):
    """
    Given a list of (x, y) pair tuples representing star positions on an image,
    return a list with the indices of the isolated stars.
    Isolated stars are stars that are at least farther than `dist_min` from any other
    star.
    Return list may be empty if no stars are isolated.
    """
    isolated = []
    kdtree = create(star_list)
    for star_index, star in enumerate(star_list):
        dists = kdtree.search_knn(star, k=2)
        if dists[1][1] > dist_min:
            isolated.append(star_index)
    return isolated
