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

        # triangulate before shifting points onto a geoid
        tri = Delaunay(points)

        # shift points to fit geoid
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
        wire_num=0
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
                        llh=coord.ecef2geodetic(r[0],r[1],r[2])
                        
                        east=coord.enu2ecef(llh[0],llh[1],0.0,1.0,0,0)
                        north=coord.enu2ecef(llh[0],llh[1],0.0,0,1.0,0)
                        
                        wires[wire_id]={"del_l":del_l,
                                        "del_e":n.dot(del_l/n.linalg.norm(del_l),east),
                                        "del_n":n.dot(del_l/n.linalg.norm(del_l),north),                                        
                                        "r":r,
                                        "llh":llh,
                                        "wire_num":wire_num,
                                        "e0":e0,
                                        "e1":e1}
                        wire_num+=1
                    else:
                        pass

                tri["wire_ids"]=wire_ids
                tris.append(tri)
            else:
                print("middle point not in tri!")

        self.llhs=llhs
        self.tris=tris
        self.points=points
        self.wires=wires
        
    def get_wire_idx(self,e0,e1):
        wire_id=self.get_wire_id(e0,e1)
        return(self.wires[wire_id]["wire_num"])
    
    def get_wire_id(self,e0,e1):
        e0p=n.min([e0,e1])
        e1p=n.max([e0,e1])
        return("%d-%d"%(e0p,e1p))
        
    def plot_currents(self,currents):
        for k in self.wires.keys():
            widx=self.wires[k]["wire_num"]
            llh=self.wires[k]["llh"]            
            I=currents[widx]
            plt.plot(llh[0],I,".")
        plt.show()

    def plot_currents2(self,currents):
        cmax=n.max(currents)
        for k in self.wires.keys():
            widx=self.wires[k]["wire_num"]
            llh=self.wires[k]["llh"]            
            I=currents[widx]
            plt.scatter(llh[1],llh[0],s=I/10,color="black")
        plt.show()
        
    def plot_currents3(self,currents):
        cmax=n.max(currents)
        c=10.0/cmax
        ax = plt.axes()
        for tri in self.tris:
            e0=tri["edges"][0]
            e1=tri["edges"][1]
            e2=tri["edges"][2]
            id0=self.get_wire_id(e0[0],e0[1])
            id1=self.get_wire_id(e1[0],e1[1])
            id2=self.get_wire_id(e2[0],e2[1])
            w0=self.wires[id0]
            w1=self.wires[id1]
            w2=self.wires[id2]            

            llh=w0["llh"]

            EI0=currents[w0["wire_num"]]
            EI1=currents[w1["wire_num"]]
            I2=currents[w2["wire_num"]]

            IE=currents[w0["wire_num"]]*w0["del_e"]+currents[w1["wire_num"]]*w1["del_e"]+currents[w2["wire_num"]]*w2["del_e"]
            IN=currents[w0["wire_num"]]*w0["del_n"]+currents[w1["wire_num"]]*w1["del_n"]+currents[w2["wire_num"]]*w2["del_n"]
            ax.arrow(llh[1],llh[0],IE*c,IN*c,head_width=0, head_length=0, color="k")
            ax.plot(llh[1],llh[0],".",color="black")
        plt.xlim([-180,180])
        plt.ylim([-90,90])        
        plt.show()
        
        
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
            ax.plot([wire["r"][0]],[wire["r"][1]],[wire["r"][2]],".",color="red")
            
        plt.show()                

if __name__ == "__main__":

    wg=wire_grid()
    wg.mean_distance()
    wg.plot()    
    wg.plot_connections()
    
    
