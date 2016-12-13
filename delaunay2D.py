"""
Simple structured Delaunay triangulation in 2D
  
Written by Jose M. Espadero ( http://github.com/jmespadero/pyDelaunay2D )
Based on code from Ayron Catteau. Published at http://github.com/ayron/delaunay

Pretend to be simple and didactic. The only requisite is numpy.
Lack of robustness checks, so do not spect to work on degenerate set of points.
"""

import numpy as np

class Triangle:
    """Simple class to store an indexed triangle with neighbour information."""

    def __init__(self, a, b, c):
        """Constructor from 3 indexes to vertex"""
        self.v = [a, b, c]
        self.neighbour = [None]*3  # Adjacent triangles
        #Avoid computing circumcenter here, so the triangle does
        #not need to know its coordinates
        self.center = None
        self.radius = None

    def __repr__(self):
        """Dump indexes as a text string"""
        return 'tri< ' + str(self.v) + ' >'
        
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
        smin= np.amin(seeds, axis=0)
        smax= np.amax(seeds, axis=0)
        #If no margin is given, compute one based on the diagonal
        if not margin:
            margin = 50*np.linalg.norm(smax- smin)
            
        # Extend the BBox coordinates by the margin. Store it.
        smin -= margin
        smax += margin
        self.margin = margin
        
        # Create a two triangles 'frame' big enough to contains all seeds
        self.coords = [smin, np.array((smax[0], smin[1])),
                       smax, np.array((smin[0], smax[1]))]

        T1 = Triangle(0, 3, 1)
        T2 = Triangle(2, 1, 3)
        T1.neighbour[0] = T2
        T2.neighbour[0] = T1
        T1.center, T1.radius = self.Circumcenter(T1)
        T2.center, T2.radius = self.Circumcenter(T2)
        self.triangles = [T1, T2]        

        #Insert all seeds one by one
        for s in seeds:
            self.AddPoint(s)
            
    def Circumcenter(self, tri):
        """Compute Circumcenter and circumradius of a triangle.
        Uses an extension of the method described here:
        http://www.ics.uci.edu/~eppstein/junkyard/circumcenter.html
        """
        pts = np.asarray([self.coords[v] for v in tri.v])              
        A = np.bmat( [[ 2*np.dot(pts,pts.T), np.ones((3,1)) ],
                   [  np.ones((1,3)) ,  np.zeros((1,1))  ]] )

        b = np.hstack((np.sum(pts * pts, axis=1), np.ones((1))))
        x = np.linalg.solve(A,b)
        bary_coords = x[:-1]  
        
        center = np.dot(bary_coords,pts)    
        radius = np.linalg.norm(pts[0] - center)
        #print("center:", center, "radius:", radius)
        return (center,radius)    

    def AddPoint(self, p):       
        """Add a new point to the current triangulation.
        """
        p = np.asarray(p)        
        idx = len(self.coords)
        #print("coords[", idx,"] ->",p)
        self.coords.append(p)
        
        # Search for the triangle(s) too near of p
        bad_triangles = []
        for T in self.triangles:            
            if np.linalg.norm(T.center - p) <= T.radius:
                bad_triangles.append(T)

        # Find the convex hull of the bad triangles.
        # Expressed as a list of edges (point pairs) in ccw order
        boundary = self.Boundary(bad_triangles)

        #Remove triangles too near of point p
        for T in bad_triangles:
            self.triangles.remove(T)

        # Retriangle the hole left by those triangles
        new_triangles = []
        for edge in boundary:
            #Create a new Triangle using point p and the edge
            T = Triangle(idx, edge[0], edge[1])
            T.center, T.radius = self.Circumcenter(T)

            #Set triangle stored at edge[2] as neighbour of T
            T.neighbour[0] = edge[2]
            
            #Set T as neighbour of triangle stored at edge[2]
            if T.neighbour[0]:
                #search the reversed edge in T.neighbour[0].v
                tmp_v = T.neighbour[0].v
                for i in range(3):
                    if edge[1] == tmp_v[i] and edge[0] == tmp_v[(i+1)%3]:
                        T.neighbour[0].neighbour[(i+2)%3] = T

            #Add triangle to a temporal list
            new_triangles.append(T)

        # Link the new triangles
        N = len(new_triangles)
        for i, T in enumerate(new_triangles):
            T.neighbour[2] = new_triangles[(i-1) % N]   # back
            T.neighbour[1] = new_triangles[(i+1) % N]   # forward
   
        #Add triangles to our triangulation
        self.triangles.extend(new_triangles)

      
    def Boundary(self, bad_triangles):
        """Build the CCW boundary of a set of triangles.
        """

        boundary = []

        #Backtracking on each triangle searching edges in the boundary
        T = bad_triangles[0]
        edge = 0
        while True:                    
            if T.neighbour[edge] in bad_triangles:
                #Move to neighbour triangle
                last = T
                T = T.neighbour[edge]
                edge = (T.neighbour.index(last) + 1) % 3 

            else:   # Found an edge that is on the boundary
                # Add edge to boundary list
                boundary.append((T.v[(edge+1)%3], T.v[(edge+2)%3], T.neighbour[edge]))
                #Move to next edge in this triangle
                edge = (edge + 1) % 3

                #Check if boundary is a closed loop
                if boundary[0][0] == boundary[-1][1]:
                    break

        return boundary

    def exportDT(self):
        """Export the current Delaunay Triangulation.
        """
        #Filter out coordinates in the extended BBox
        xs = [p[0] for p in self.coords[4:] ]
        ys = [p[1] for p in self.coords[4:] ]        
        
        #Filter out triangles with any vertex in the extended BBox
        ETS = [t.v for t in self.triangles ]
        tris = [(a-4,b-4,c-4) for (a,b,c) in ETS if a > 3 and b > 3 and c > 3]
        return xs, ys, tris
        
    def exportExtendedDT(self):
        """Export the Extended Delaunay Triangulation (include the frame vertex).
        """
        xs = [p[0] for p in self.coords ]
        ys = [p[1] for p in self.coords ]        
        tris = [t.v for t in self.triangles ]
        return xs, ys, tris

