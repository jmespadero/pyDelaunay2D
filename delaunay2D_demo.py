#!/usr/bin/env python3
# -*- coding: ascii -*-
"""
Simple delaunay2D demo with mathplotlib
Written by Jose M. Espadero < http://github.com/jmespadero/pyDelaunay2D >
"""
import numpy as np
from delaunay2D import Delaunay2D

if __name__ == '__main__':

    ###########################################################
    # Generate 'numSeeds' random seeds in a square of size 'radius'
    numSeeds = 24
    radius = 100
    seeds = radius * np.random.random((numSeeds, 2))
    print("seeds:\n", seeds)
    print("BBox Min:", np.amin(seeds, axis=0),
          "Bbox Max: ", np.amax(seeds, axis=0))

    """
    Compute our Delaunay triangulation of seeds.
    """
    # It is recommended to build a frame taylored for our data
    # dt = D.Delaunay2D() # Default frame
    dt = Delaunay2D(center=np.mean(seeds, axis=0), radius = 50* radius)

    # Insert all seeds one by one
    for s in seeds:
        dt.AddPoint(s)
    
    # Extract our result (without extended BBox coordinates)
    dt_x, dt_y, dt_tris = dt.exportDT()
    print("dt_tris:", len(dt_tris), "triangles")
    print("dt_tris:\n", dt_tris)
    # print triangles with neighbours
    # for t in dt.triangles:        
    #     print(t, dt.triangles[t])
    
    """
    Show how to plot triangular grids.
    """
    import matplotlib.pyplot as plt
    import matplotlib.tri
    import matplotlib.collections

    # Create a plot with matplotlib.pyplot
    fig, ax = plt.subplots()
    ax.margins(0.1)
    ax.set_aspect('equal')

    # Plot our Delaunay triangulation (plot in blue)
    ax.triplot(matplotlib.tri.Triangulation(dt_x, dt_y, dt_tris), 'bo--')

    # DEBUG: Use matplotlib to create a Delaunay triangulation (plot in green)
    # DEBUG: It should be equal to our result in dt_tris (plot in blue)
    # DEBUG: If boundary is diferent, try to increase the value of your margin
    # ax.triplot(matplotlib.tri.Triangulation(dt_x, dt_y), 'g--')

    # Plot the circumcircles (circles in black)
    # for c in dt.exportCircles():
    #   ax.add_artist(plt.Circle(c[0], c[1], color='k', fill=False, ls='dotted'))

    # Plot voronoi diagram edges (in red)
    ve = dt.exportVoronoiEdges()
    ax.add_collection(matplotlib.collections.LineCollection(ve, colors='r'))
    
    # DEBUG: plot the extended triangulation (plot in red)
    # edt_x, edt_y, edt_tris = dt.exportExtendedDT()
    # print("Extended dt_tris:", len(edt_tris), "triangles")
    # ax.triplot(matplotlib.tri.Triangulation(edt_x, edt_y, edt_tris), 'ro-.')

    # Dump plot to file
    # plt.savefig('output-delaunay2D.png', dpi=96)
    # plt.savefig('output-delaunay2D.svg', dpi=96)

    plt.show()

    """
    Demo of a step-by-step triangulation plot
    """

    """
    # Build another DT frame
    dt2 = Delaunay2D(center=np.mean(seeds, axis=0), radius = 100* radius)
    for i,s in enumerate(seeds):
        dt2.AddPoint(s)
        if i > 1:
            #plt.gcf().clear()
            fig, ax = plt.subplots()
            ax.margins(0.1)
            ax.set_aspect('equal')
            dt_x, dt_y, dt_tris = dt2.exportDT()
            ax.triplot(matplotlib.tri.Triangulation(dt_x, dt_y, dt_tris))
            plt.show()
    """
