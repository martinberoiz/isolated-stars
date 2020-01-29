# -*- coding: utf-8 -*-
"""
Code Challenge: Isolated Stars.
(c) CGWA Time Domain Astronomy Group
"""
from scipy.spatial import KDTree

def isolated(star_list, dist_min):
    """
    Given a list of (x, y) pair tuples representing star positions on an image,
    return a list with the indices of the isolated stars.
    Isolated stars are stars that are at least farther than `dist_min` from any other
    star.
    Return list may be empty if no stars are isolated.
    """
    isolated = []
    kdtree = KDTree(star_list)
    for star_index, star in enumerate(star_list):
        dists, close_stars = kdtree.query(star, k=2)
        if dists[1] > dist_min:
            isolated.append(star_index)
    return isolated
