#!/usr/bin/env python

import numpy as n

import sphere as sg
import supermag_read as sr
import coord
import scipy.constants as sc
import matplotlib.pyplot as plt

def create_theory_matrix(meas,h_e=100e3,R_e=6378e3,min_lat=60,max_lat=80,min_lon=3,max_lon=30,n_lat_wires=40,n_seg=50):
    """
    Given measurements, determine the current pattern at h_e km
    """
    lats=n.linspace(min_lat,max_lat,num=n_lat_wires+1)
    lons_seg=n.linspace(min_lon,max_lon,num=n_seg+1)

    wlens=n.zeros([n_lat_wires,n_seg])
    rl=[]
    dls=[]
    dAs=[]
    for wi in range(n_lat_wires):
        r=[]
        dl=[]
        dA=[]
        for si in range(n_seg):
            
            l0=coord.geodetic2ecef(lats[wi],lons_seg[si],h_e)
            l0u=coord.geodetic2ecef(lats[wi+1],lons_seg[si],h_e)
            
            l1=coord.geodetic2ecef(lats[wi],lons_seg[si+1],h_e)
            l1u=coord.geodetic2ecef(lats[wi+1],lons_seg[si+1],h_e)

            # area in this pixel
            A=n.linalg.norm(l0u-l0)*n.linalg.norm(l1u-l1)
            
            del_l=(l1-l0)
            # unit vector
            del_l=del_l/n.linalg.norm(del_l)
            
            dl.append(del_l)
            r.append(0.5*(l0+l1))
            dA.append(A)
            
        rl.append(r)
        dls.append(dl)
        dAs.append(dA)
    lats=lats[0:(len(lats)-1)]
    n_stat=len(meas["glon"])
    n_meas=3*n_stat

    n_par=n_lat_wires
    print("n_meas %d n_par %d"%(n_meas,n_par))
    
    # theory matrix
#    A=n.zeros([n_meas,n_par])
    # measurements
 #   m=n.zeros(n_meas)
    
    const=sc.mu_0/4.0/n.pi
    A={}
    for mi in range(n_stat):
        p=coord.geodetic2ecef(meas["glat"][mi],meas["glon"][mi],0)
        
        up=coord.enu2ecef(meas["glat"][mi],meas["glon"][mi],0,0,0,1)
        Ax=n.zeros(n_lat_wires)
        Ay=n.zeros(n_lat_wires)
        Az=n.zeros(n_lat_wires)        
        for li in range(n_lat_wires):
            for si in range(n_seg):
                # for each wire segment, get the Biot-Savart law contribution
                # is ionosphere above the horizon?
                rp = rl[li][si]-p
                rpn=n.linalg.norm(rp)
                rp0 = rp/rpn
                zen_angle=180.0*n.arccos(n.dot(rp0,up))/n.pi
                if zen_angle < 85.0:
                    
                    print("station %s lat %1.2f lon %1.2f grid lat %1.2f lon %1.2f"%(meas["stat"][mi],
                                                                                     meas["glat"][mi],
                                                                                     meas["glon"][mi],
                                                                                     lats[li],
                                                                                     lons_seg[si]))
                    # biot-savart law (current density per unit surface area)
                    cp=dAs[li][si]*const*n.cross(dls[li][si],rp)/rpn**3.0
                    # x
                    Ax[li]+=cp[0]
                    # y
                    Ay[li]+=cp[1]                
                    # z
                    Az[li]+=cp[2]
        if n.sum(Ax) != 0:
            stat=meas["stat"][mi]
            A[stat]={"x":Ax,"y":Ay,"z":Az}
    n_stations=len(A.keys())
    print("n_stat %d"%(n_stations))
    return(A,lats,lons_seg)


def populate_theory_matrix(A,meas,reg=1e-4):
    stats=A.keys()

    Alist=[]
    mlist=[]
    n_meas=0
    for s in stats:
        midx=n.where(meas["stat"] == s)[0]
        for mi in midx:
            # convert the n,e,u meas to ecef in Tesla                             E                N           U
            B_meas=coord.enu2ecef(meas["glat"][mi],meas["glon"][mi],0,meas["neu"][mi,1],meas["neu"][mi,0],meas["neu"][mi,2])*1e-9
            Alist.append(A[s]["x"])
            Alist.append(A[s]["y"])
            Alist.append(A[s]["z"])
            mlist.append(B_meas[0])
            mlist.append(B_meas[1])
            mlist.append(B_meas[2])
            n_meas+=3
    n_par=len(Alist[0])
    for ri in range(n_par):
        z=n.zeros(n_par)
        z[ri]=reg
        Alist.append(z)
        mlist.append(0.0)
    AM=n.array(Alist)
    m=n.array(mlist)
        
    return(AM,m,n_meas)

    
d=sr.supermag_data(fname="test_bl_l_30.csv")
b=d.get_bounds()

meas=d.get_meas(b[0],b[0]+60)
print(n.mean(meas["times"]))

n_lat_wires=20
Ak,lats,lons=create_theory_matrix(meas,n_lat_wires=n_lat_wires)

dt=120.0
Nt=int(24*3600*7/dt)
Ir=n.zeros([Nt,n_lat_wires])
I=n.zeros([Nt,n_lat_wires])
times=[]
stds=[]
for mi in range(Nt):
    times.append(b[0]+mi*dt)
    meas=d.get_meas(b[0]+mi*dt,b[0]+60+mi*dt)
    A,m,n_meas=populate_theory_matrix(Ak,meas)
    
    
    xhat_ls=n.linalg.lstsq(A,m)[0]
    model=n.dot(A,xhat_ls)
#    plt.plot((m[0:n_meas]-model[0:n_meas])*1e9)
    s=n.std((m[0:n_meas]-model[0:n_meas])*1e9)
    stds.append(s)
  #  plt.show()



    if False:
        Sigma_P=mvar*n.linalg.inv(n.dot(n.transpose(A),A))        
        plt.plot(xhat_ls)
        plt.plot(xhat_ls+n.sqrt(n.diag(Sigma_P)),color="grey")
        plt.plot(xhat_ls-n.sqrt(n.diag(Sigma_P)),color="grey")
        plt.show()
    
#    plt.plot(m)
#    plt.plot(m-model)
#    plt.show()
    
    I[mi,:]=xhat_ls

    u,s,vh=n.linalg.svd(A)

    sinv=n.zeros([A.shape[1],u.shape[0]])
 #   plt.plot(s)
#    plt.show()
    # s/(s**2
#    bidx=n.where(s < n.max(s)/2.0)[0]
    sinv0=s/(s**2+0.25**2*n.max(s)**2.0)
#    sinv0[bidx]=0.0

    #plt.semilogy(sinv0)
    #plt.semilogy(1.0/s)
    #plt.show()
    print(sinv.shape)
    #else:
    #   sinv=n.zeros([len(s),u.shape[0]])
    
    for i in range(n.min([sinv.shape[0],sinv.shape[1]])):
        sinv[i,i]=sinv0[i]
        
    xhat=n.dot(n.transpose(vh),n.dot(sinv,n.dot(n.transpose(u),m)))
    Ir[mi,:]=xhat
plt.plot(stds)
plt.show()
plt.pcolormesh((times-times[0])/3600.0,lats,n.transpose(I)*1e3*1e3)
plt.colorbar()
plt.show()
plt.pcolormesh((times-times[0])/3600.0,lats,n.transpose(Ir)*1e3*1e3)
plt.colorbar()
plt.show()




