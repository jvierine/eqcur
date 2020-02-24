#!/usr/bin/env python

import numpy as n

import sphere as sg
import supermag_read as sr
import coord
import scipy.constants as sc
import matplotlib.pyplot as plt

def create_theory_matrix(wg,meas,h_e=100e3,R_e=6378e3,alpha=1e-11,tik=1e-13):
    """
    Given measurements, determine the current pattern at h_e km
    """
    n_stat=len(meas["glon"])
    n_meas=3*n_stat

    n_par=len(wg.wires)
    print("n_meas %d n_par %d"%(n_meas,n_par))    
    # how many points can we regularize using current continuity
    n_cur_reg=wg.points.shape[0]-1
    n_tik=n_par
    n_reg=n_cur_reg#+n_tik
    
    # theory matrix
    A=n.zeros([n_meas+n_reg,n_par])
    # measurements
    m=n.zeros(n_meas+n_reg)
    
    wire_ids=wg.wires.keys()
    
    const=sc.mu_0/4.0/n.pi
    
    for mi in range(n_stat):
        p=coord.geodetic2ecef(meas["glat"][mi],meas["glon"][mi],0)
        # convert the n,e,u meas to ecef in Tesla                             E                N           U
        B_meas=coord.enu2ecef(meas["glat"][mi],meas["glon"][mi],0,meas["neu"][mi,1],meas["neu"][mi,0],meas["neu"][mi,2])*1e-9

        up=coord.enu2ecef(meas["glat"][mi],meas["glon"][mi],0,0,0,1)

        for li in range(n_par):
            k=wire_ids[li]
            # for each wire, get the Biot-Savart law contribution
            # is ionosphere above the horizon?
            rp = wg.wires[k]["r"]-p
            rpn=n.linalg.norm(rp)
            rp0 = rp/rpn
            zen_angle=180.0*n.arccos(n.dot(rp0,up))/n.pi
            if zen_angle < 90.0:
#                print("station lat %1.2f lon %1.2f grid lat %1.2f lon %1.2f"%(meas["glat"][mi],
 #                                                                             meas["glon"][mi],
  #                                                                            wg.wires[k]["llh"][0],
   #                                                                           wg.wires[k]["llh"][1]))
                # biot-savart law
                cp=const*n.cross(wg.wires[k]["del_l"],rp)/rpn**3.0
                wire_idx=wg.wires[k]["wire_num"]
                # x
                A[mi*3+0,wire_idx]=cp[0]
                # y
                A[mi*3+1,wire_idx]=cp[1]                
                # z
                A[mi*3+2,wire_idx]=cp[2]
                # measurement
                m[mi*3+0]=B_meas[0]
                m[mi*3+1]=B_meas[1]
                m[mi*3+2]=B_meas[2]
                
        # current continuity
        nodes=wg.connections.keys()
        for ri in range(n_cur_reg):
            e0=nodes[ri]
            conns=wg.connections[e0]
            for e1 in conns:
                wire_idx=wg.get_wire_idx(e0,e1)
                # if index is reversed, current negative
                if e1 > e0:
                    A[ri+n_stat*3,wire_idx]=alpha
                else:
                    A[ri+n_stat*3,wire_idx]=-alpha
                    
#        for ri in range(n_tik):
 #           A[ri+n_stat*3+n_cur_reg,ri]=tik

    return(A,m)
    print("n_meas %d n_reg %d n_par %d"%(n_meas,n_reg,n_par))
    
# create a wiregrid model
wg=sg.wire_grid(N=1000)
wg.plot()
d=sr.supermag_data()
b=d.get_bounds()
h=19.27
meas=d.get_meas(b[0]+h*3600,b[0]+60+h*3600)

A,m=create_theory_matrix(wg,meas)
u,s,vh=n.linalg.svd(A)
print(A.shape)
print(len(m))
print(u.shape)
print(vh.shape)
print(len(s))
#sinv=1.0/s
#if len(s) > u.shape[0]:


sinv=n.zeros([A.shape[1],u.shape[0]])
# s/(s**2
bidx=n.where(s < n.max(s)/2.0)[0]
sinv0=1.0/s
sinv0[bidx]=0.0

plt.semilogy(sinv0)
plt.semilogy(1.0/s)
plt.show()
print(sinv.shape)
#else:
 #   sinv=n.zeros([len(s),u.shape[0]])
    
for i in range(n.min([sinv.shape[0],sinv.shape[1]])):
    sinv[i,i]=sinv0[i]
    
xhat=n.dot(n.transpose(vh),n.dot(sinv,n.dot(n.transpose(u),m)))


#plt.pcolormesh(1e12*A)
#plt.colorbar()
#plt.show()

#xhat=n.linalg.lstsq(1e12*A,1e12*m)[0]
model=n.dot(A,xhat)
n_meas=len(meas["glon"])*3

plt.subplot(121)
plt.plot(m-model)
plt.subplot(122)
plt.plot(m)
plt.show()

wg.plot_currents3(xhat)

