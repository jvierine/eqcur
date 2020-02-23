#!/usr/bin/env python

import numpy as n

import sphere as sg
import supermag_read as sr
import coord

def create_theory_matrix(wg,meas,h_e=100e3,R_e=6378e3):
    """
    Given measurements, determine the current pattern at h_e km
    """
    n_meas=len(meas["glon"])
    n_par=len(wg.tris)
    print("n_meas %d n_par %d"%(n_meas,n_par))
#    A=n.zeros([n_meas*3,

# create a wiregrid model
wg=sg.wire_grid(N=500)
wg.plot()

d=sr.supermag_data()
b=d.get_bounds()
meas=d.get_meas(b[0],b[0]+60)

A=create_theory_matrix(wg,meas)
