#!/usr/bin/env python

import numpy as n
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial import Delaunay
import supermag_stations as ss

import math, random

def fibonacci_sphere(samples=1,randomize=False):
    """
    Somewhat uniform grid on a sphere
    """
    rnd = 1.
    if randomize:
        rnd = random.random() * samples

    points = []
    points.append([0,0,0])
    offset = 2./samples
    increment = math.pi * (3. - math.sqrt(5.));

    for i in range(samples):
        y = ((i * offset) - 1) + (offset / 2);
        r = math.sqrt(1 - pow(y,2))

        phi = ((i + rnd) % samples) * increment

        x = math.cos(phi) * r
        z = math.sin(phi) * r

        points.append([x,y,z])

    points=n.array(points)
    return(points)

class wire_grid():
    """
    Let's do some geometry here. 
    We assume that the ionosphere is a grid of current
    carrying wires. 
    """
    def __init__(self,N=1000,R_e=6378e3,plot=False):
        # get uniformly distributed points on a sphere
        # these will act ast grid points
        points=fibonacci_sphere(samples=N)*R_e
        tri = Delaunay(points)
        # triangulate the grid, so that we know
        # the neighbouring points
        # the list "tri" contains the wire mesh

        # TBD: make sure that there are no duplicated wires.
            
        s=n.copy(tri.simplices)
        tris = []
        
        for i in range(s.shape[0]):
            pidx=s[i,:]
            
            if 0 in pidx:
                pidx=n.setdiff1d(pidx,[0])
                
                tri={}
                
                tri["edges"]=[[pidx[0],pidx[1]],[pidx[0],pidx[2]],[pidx[1],pidx[2]]]
                
                tri["del_l"]=[ points[pidx[0],:]-points[pidx[1],:],
                               points[pidx[0],:]-points[pidx[2],:],
                               points[pidx[1],:]-points[pidx[2],:] ]
                
                tri["r"]=[ points[pidx[0],:]+0.5*(points[pidx[0],:]-points[pidx[1],:]),
                           points[pidx[0],:]+0.5*(points[pidx[0],:]-points[pidx[2],:]),
                           points[pidx[1],:]+0.5*(points[pidx[1],:]-points[pidx[2],:]) ]
                tris.append(tri)
                
            else:
                print("middle point not in tri!")
        self.tris=tris
        self.points=points
        
    def plot(self):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        
        for tri in self.tris:
            for e in tri["edges"]:
                ax.plot([self.points[e[0],0],self.points[e[1],0]],
                        [self.points[e[0],1],self.points[e[1],1]],
                        [self.points[e[0],2],self.points[e[1],2]],color="grey")

        plt.show()                
        
        
                



if __name__ == "__main__":

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    wg=wire_grid()
    s=ss.supermag_stations()
    stats=s.props.keys()

    for tri in wg.tris:
        for e in tri["edges"]:
            ax.plot([wg.points[e[0],0],wg.points[e[1],0]],
                    [wg.points[e[0],1],wg.points[e[1],1]],
                    [wg.points[e[0],2],wg.points[e[1],2]],color="grey")
            
    
    for sn in stats:

        p=s.prop(sn)
        print("%s lat %1.2f lon %1.2f"%(sn,p["glat"],p["glon"]))

        phi=n.pi*p["glat"]/180.0
        theta=n.pi*p["glon"]/180.0

        x=n.cos(phi)*n.cos(theta)
        y=n.cos(phi)*n.sin(theta)
        z=n.sin(phi)
        ax.plot([x],[y],[z],".",color="red")
        
    plt.show()

    for sn in stats:
        p=s.prop(sn)
        print("%s lat %1.2f lon %1.2f"%(sn,p["glat"],p["glon"]))
        plt.plot(p["glon"],p["glat"],".")
        
    plt.show()
    
