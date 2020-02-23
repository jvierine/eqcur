#!/usr/bin/env python

import numpy as n
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial import Delaunay
import supermag_stations as ss

import coord

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
    def __init__(self,N=1000,R_e=6378e3,h_e=100e3,plot=False):
        self.N=N
        self.R_e=R_e
        # get uniformly distributed points on a sphere
        # these will act ast grid points
        points=fibonacci_sphere(samples=N)*(R_e+h_e)
        
        tri = Delaunay(points)
        
        # fit geoid
        if True:
            points2=n.copy(points)
            llhs=[]
            for p in range(points.shape[0]):
                if p == 0:
                    llh=[0,0,0]
                    points[p,:]=n.array([0,0,0])
                else:
                    llh=coord.ecef2geodetic(points[p,0],points[p,1],points[p,2])
                    points2[p,:]=coord.geodetic2ecef(llh[0],llh[1],h_e)
                llhs.append(llh)
            llhs=n.array(llhs)
            points=points2
        
        # triangulate the grid, so that we know
        # the neighbouring points
        # the list "tri" contains the wire mesh

        # TBD: make sure that there are no duplicated wires!
            
        s=n.copy(tri.simplices)

        # these are the triangles
        tris = []
        
        # these are unique current carrying elements of the mesh
        wires={}
        
        # these are the connections between node points, in order to allow us to regularize current continuity.
        self.connections={}
        
        for i in range(s.shape[0]):
            pidx=s[i,:]
            
            if 0 in pidx:
                pidx=n.setdiff1d(pidx,[0])
                
                tri={}
                
                tri["edges"]=[[pidx[0],pidx[1]],[pidx[0],pidx[2]],[pidx[1],pidx[2]]]

                id_pairs=[[0,1],[0,2],[1,2]]

                wire_ids=[]
                for id_pair in id_pairs:
                    
                    e0=n.min([pidx[id_pair[0]],pidx[id_pair[1]]])
                    e1=n.max([pidx[id_pair[0]],pidx[id_pair[1]]])

                    self.add_connection(e0,e1)
                    
                    wire_id = "%d-%d"%(e0,e1)
                    wire_ids.append(wire_id)
                    if wire_id not in wires.keys():
                        del_l=points[e0,:]-points[e1,:]
                        r=points[e0,:]+0.5*del_l
                        wires[wire_id]={"del_l":del_l,
                                        "r":r,
                                        "e0":e0,
                                        "e1":e1}
                    else:
                        pass

                tri["wire_ids"]=wire_ids
                tris.append(tri)
            else:
                print("middle point not in tri!")
                
        self.tris=tris
        self.points=points
        self.wires=wires
        print(len(self.points))        
        print(len(self.wires.keys()))

    def add_connection(self,e0,e1):
        """
        In order to treat the grid as a circuit, we need to 
        know how points are connected.
        when connecting from smaller idx to larger idx, the current is positive, otherwise negative 
        """
        if e0 in self.connections.keys():
            if e1 not in self.connections[e0]:
                self.connections[e0].append(e1)
        else:
            self.connections[e0]=[e1]
        # reverse connection
        if e1 in self.connections.keys():
            if e0 not in self.connections[e1]:
                self.connections[e1].append(e0)
        else:
            self.connections[e1]=[e0]
    
    def mean_distance(self):
        """
        calculate the stats on distance between nodes in grid
        """
        dists=[]
        for tri in self.tris:
            e=tri["edges"]
            dists.append(n.linalg.norm(self.points[e[0][0],:]-self.points[e[0][1],:]))
            dists.append(n.linalg.norm(self.points[e[1][0],:]-self.points[e[1][1],:]))
            dists.append(n.linalg.norm(self.points[e[2][0],:]-self.points[e[2][1],:]))
        print("distance stats")
        
        sanity_check=n.sqrt((4*n.pi*self.R_e**2.0)/self.N)
        print("mean length when divided into squares %1.2f km"%(sanity_check/1e3))
        print(n.mean(dists)/1e3)
        print(n.min(dists)/1e3)
        print(n.max(dists)/1e3)
            
    def plot_connections(self):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        
        for k in self.connections.keys():
            
            print("%d connected to"%(k))
            print(self.connections[k])
            for c in self.connections[k]:
                ax.plot([self.points[k,0],self.points[c,0]],
                        [self.points[k,1],self.points[c,1]],
                        [self.points[k,2],self.points[c,2]])
        ax.set_xlim([-7000e3,7000e3])
        ax.set_ylim([-7000e3,7000e3])
        ax.set_zlim([-7000e3,7000e3])            
        plt.show()
        
    def plot(self):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        for wk in self.wires.keys():
            wire=self.wires[wk]
            e0=wire["e0"]
            e1=wire["e1"]            
            ax.plot([self.points[e0,0],self.points[e1,0]],
                    [self.points[e0,1],self.points[e1,1]],
                    [self.points[e0,2],self.points[e1,2]],color="grey")
#            print(coord.ecef2geodetic(wire["r"][0],wire["r"][1],wire["r"][2]))
            ax.plot([wire["r"][0]],[wire["r"][1]],[wire["r"][2]],".",color="red")
            
        plt.show()                

if __name__ == "__main__":

    wg=wire_grid()
    wg.mean_distance()
    wg.plot()    
    wg.plot_connections()
    
    
