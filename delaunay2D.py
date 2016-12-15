# -*- coding: ascii -*-
"""
Simple structured Delaunay triangulation in 2D with Bowyer-Watson algorithm.

Written by Jose M. Espadero ( http://github.com/jmespadero/pyDelaunay2D )
Based on code from Ayron Catteau. Published at http://github.com/ayron/delaunay

Just pretend to be simple and didactic. The only requisite is numpy.
Robust checks disabled by default. May not work in degenerate set of points.
"""

import numpy as np
from math import sqrt


class Delaunay2D:
    """
    Class to compute a Delaunay triangulation in 2D
    ref: http://en.wikipedia.org/wiki/Bowyer%E2%80%93Watson_algorithm
    ref: http://www.geom.uiuc.edu/~samuelp/del_project.html
    """

    def __init__(self, seeds, margin=None):
        """Create a new set of triangles from a 2D vertex array
        seeds  -- Array of 2D coordinates for the vertex
        margin -- Optional parameter to set the distance of extended BBox
        """

        seeds = np.asarray(seeds)

        # Compute BBox
        smin = np.amin(seeds, axis=0)
        smax = np.amax(seeds, axis=0)
        # If no margin is given, compute one based on the diagonal
        if not margin:
            margin = 50 * np.linalg.norm(smax - smin)

        # Extend the BBox coordinates by the margin. Store it.
        smin -= margin
        smax += margin
        self.margin = margin

        # Create a two triangles 'frame' big enough to contains all seeds
        self.coords = [smin, np.array((smax[0], smin[1])),
                       smax, np.array((smin[0], smax[1]))]

        # Create two triangles for the frame
        T1 = (0, 3, 1)
        T2 = (2, 1, 3)
        self.triangles = {T1: [T2, None, None], T2: [T1, None, None]}

        # Compute circumcenters and circumradius for each triangle
        self.circles = {}
        for t in self.triangles:
            self.circles[t] = self.Circumcenter(t)

        # Insert all seeds one by one
        for s in seeds:
            self.AddPoint(s)

    def Circumcenter(self, tri):
        """Compute Circumcenter and circumradius of a triangle in 2D.
        Uses an extension of the method described here:
        http://www.ics.uci.edu/~eppstein/junkyard/circumcenter.html
        """
        pts = np.asarray([self.coords[v] for v in tri])
        pts2 = np.dot(pts, pts.T)
        A = np.bmat([[2 * pts2, [[1],
                                 [1],
                                 [1]]],
                      [[[1, 1, 1, 0]]]])

        b = np.hstack((np.sum(pts * pts, axis=1), [1]))
        x = np.linalg.solve(A, b)
        bary_coords = x[:-1]
        center = np.dot(bary_coords, pts)

        # radius = np.linalg.norm(pts[0] - center) # euclidean distance
        radius = np.sum(np.square(pts[0] - center))  # squared distance
        return (center, radius)

    def inCircleFast(self, tri, p):
        """Check if point p is inside of precomputed circumcircle of tri.
        """
        circle = self.circles[tri]
        return np.sum(np.square(circle[0] - p)) <= circle[1]

    def inCircleRobust(self, tri, p):
        """Check if point p is inside of circumcircle around the triangle tri.
        This is a robust predicate, slower than compare distance to centers
        ref: http://www.cs.cmu.edu/~quake/robust.html
        """
        m1 = np.asarray([self.coords[v] - p for v in tri])
        m2 = np.sum(np.square(m1), axis=1).reshape((3, 1))
        m = np.hstack((m1, m2))    # The 3x3 matrix to check
        return (np.linalg.det(m) <= 0)

    def AddPoint(self, p):
        """Add a new point to the current triangulation.
        """
        p = np.asarray(p)
        idx = len(self.coords)
        # print("coords[", idx,"] ->",p)
        self.coords.append(p)

        # Search the triangle(s) whose circumcircle contains p
        bad_triangles = []
        for T in self.triangles:
            # Choose one method: inCircleRobust(T, p) or inCircleFast(T, p)
            if self.inCircleFast(T, p):
                bad_triangles.append(T)

        # Find the CCW boundary (star shape) of the bad triangles,
        # expressed as a list of edges (point pairs) and the opposite
        # triangle to each edge.
        boundary = []
        # Choose a "random" triangle and edge
        T = bad_triangles[0]
        edge = 0
        # get the opposite triangle of this edge
        while True:
            # Check if edge of triangle T is on the boundary...
            # if opposite triangle of this edge is external to the list
            tri_op = self.triangles[T][edge]
            if tri_op not in bad_triangles:
                # Insert edge and external triangle into boundary list
                boundary.append((T[(edge+1) % 3], T[(edge-1) % 3], tri_op))

                # Move to next CCW edge in this triangle
                edge = (edge + 1) % 3

                # Check if boundary is a closed loop
                if boundary[0][0] == boundary[-1][1]:
                    break
            else:
                # Move to next CCW edge in opposite triangle
                edge = (self.triangles[tri_op].index(T) + 1) % 3
                T = tri_op

        # Remove triangles too near of point p of our solution
        for T in bad_triangles:
            del self.triangles[T]
            del self.circles[T]

        # Retriangle the hole left by bad_triangles
        new_triangles = []
        for (e0, e1, tri_op) in boundary:
            # Create a new triangle using point p and edge extremes
            T = (idx, e0, e1)

            # Store circumcenter and circumradius of the triangle
            self.circles[T] = self.Circumcenter(T)

            # Set opposite triangle of the edge as neighbour of T
            self.triangles[T] = [tri_op, None, None]

            # Try to set T as neighbour of the opposite triangle
            if tri_op:
                # search the neighbour of tri_op that use edge (e1, e0)
                for i, neigh in enumerate(self.triangles[tri_op]):
                    if neigh:
                        if e1 in neigh and e0 in neigh:
                            # change link to use our new triangle
                            self.triangles[tri_op][i] = T

            # Add triangle to a temporal list
            new_triangles.append(T)

        # Link the new triangles each another
        N = len(new_triangles)
        for i, T in enumerate(new_triangles):
            self.triangles[T][1] = new_triangles[(i+1) % N]   # next
            self.triangles[T][2] = new_triangles[(i-1) % N]   # previous

    def exportDT(self):
        """Export the current Delaunay Triangulation.
        """
        # Filter out coordinates in the extended BBox
        xs = [p[0] for p in self.coords[4:]]
        ys = [p[1] for p in self.coords[4:]]

        # Filter out triangles with any vertex in the extended BBox
        tris = [(a-4, b-4, c-4)
                for (a, b, c) in self.triangles if a > 3 and b > 3 and c > 3]
        return xs, ys, tris

    def exportExtendedDT(self):
        """Export the Extended Delaunay Triangulation (with the frame vertex).
        """
        xs = [p[0] for p in self.coords]
        ys = [p[1] for p in self.coords]
        tris = [t for t in self.triangles]
        return xs, ys, tris

    def exportCircles(self):
        """Export the circumcircles
        """
        # Remember to compute circumcircles if not done before
        # for t in self.triangles:
        #     self.circles[t] = self.Circumcenter(t)

        # Filter out triangles with any vertex in the extended BBox
        # Do sqrt of radius before of return
        return [(self.circles[(a, b, c)][0], sqrt(self.circles[(a, b, c)][1]))
                for (a, b, c) in self.triangles if a > 3 and b > 3 and c > 3]

    def exportVoronoiEdges(self):
        """Export the edges of Voronoi diagram as unstructured data.
           May contain duplicates edges.
        """
        # Remember to compute circumcircles if not done before
        # for t in self.triangles:
        #     self.circles[t] = self.Circumcenter(t)
        vor_edges = []
        for (a, b, c) in self.triangles:
            if a > 3 and b > 3 and c > 3:
                t = (a, b, c)
                for n in self.triangles[t]:
                    vor_edges.append([self.circles[t][0], self.circles[n][0]])
        return vor_edges

    def exportVoronoiCoordEdges(self):
        """Export the edges of Voronoi diagram as indexed data.
           May contain duplicates
        """
        # Remember to compute circumcircles if not done before
        # for t in self.triangles:
        #     self.circles[t] = self.Circumcenter(t)
        ets = [t for t in self.triangles]
        vor_coors = [self.circles[t][0] for t in ets]
        vor_edges = []
        for tidx, (a, b, c) in enumerate(ets):
            if a > 3 and b > 3 and c > 3:
                for neigh in self.triangles[(a, b, c)]:
                    vor_edges.append((tidx, ets.index(neigh)))
        return vor_coors, vor_edges
