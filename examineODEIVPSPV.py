# -*- coding: utf-8 -*-
"""
Created on Mon Nov 24 13:31:44 2014

@author: tai
"""

import numpy as np
import pylab as pl
import pyODEIVPSPV as ode

print "#"*30
print ode.__doc__

print "#"*30
print ode.fdxsample1.__doc__

print "#"*30
print ode.fdxsample2.__doc__

print "#"*30
print ode.rk4.__doc__


# test sample state derivative functions
#px = np.array([-1,0.2,0.0,  -2,0.5,0.0,0.7,0.8], np.float64)
p = np.array([1,0.2,3.0], np.float64);
x0 = np.array([1,2], np.float64)
u0 = np.array([4], np.float64)
t0 = 0.0
"""
"""
dx = ode.fdxsample2(t0,x0,u0,p)
print dx


nu = 1
ny = 2

def input1(t,x,p): 
    return np.array([3.*np.sin(t)], np.float32)

def output1(t,x,u,p,dx): 
    return np.array([x[0]+x[1],x[0]-x[1]], np.float32)


tv = np.arange(0.0, 40.0, 0.1)
xv,yv = ode.rk4(ode.fdxsample2, input1, output1, x0, tv, p, nu, ny)

pl.figure()
pl.plot(tv,xv[:,0],tv,xv[:,1])

pl.figure()
pl.plot(tv,yv[:,0],tv,yv[:,1],
        tv,xv[:,0]+xv[:,1],'.',tv,xv[:,0]-xv[:,1])
        
pl.figure()
pl.plot(tv,yv[:,0]-(xv[:,0]+xv[:,1]),tv,yv[:,1]-(xv[:,0]-xv[:,1]))
