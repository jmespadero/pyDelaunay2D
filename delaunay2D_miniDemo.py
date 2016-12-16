#!/usr/bin/env python3
"""
Minimal delaunay2D test
See: http://github.com/jmespadero/pyDelaunay2D
"""
import numpy as np
from delaunay2D import Delaunay2D

# Create a random set of 2D points
seeds = np.random.random((10, 2))

# Create delaunay Triangulation
dt = Delaunay2D()
for s in seeds:
    dt.AddPoint(s)

# Dump triangles 
print (dt.exportTriangles())
