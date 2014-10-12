# -*- coding: utf-8 -*-
"""
Created on Sat Oct 11 17:59:12 2014

@author: hakan
"""

import pyODEIVPstage0 as ode
import numpy as np
import pylab as pl

#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
pl.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
pl.rc('text', usetex=True)

plotsOn = True


print "-"*40
print "-"*40
print "-"*40

print ode.__doc__

print "-"*40
print "-"*40
print ode.eulerintegration.__doc__
print "-"*40
print "-"*40
print ode.rk4.__doc__
print "-"*40
print "-"*40
print ode.statederivative.__doc__
print "-"*40
print "-"*40
print ode.fdxsample1.__doc__
print "-"*40
print "-"*40
print ode.fdxsample2.__doc__
print "-"*40
print "-"*40

# test sample state derivative functions
#px = np.array([-1,0.2,0.0,  -2,0.5,0.0,0.7,0.8], np.float64)
px = np.array([1,0.2,3.0], np.float64);
x0 = np.array([1,2], np.float64)
u0 = np.array([4], np.float64)
t0 = 0.0
"""
"""
dx = ode.fdxsample2(t0,x0,u0,px)
print dx

# test integrators with callbacks defined in Python
t1 = 30.0
dt = 0.1
pu = np.array([1], np.float64)
py = np.array([1], np.float64)

def fu(t,x,pu):
    delay = pu[0]
    u = [10*np.sin(t-delay)]
    #print "u is ",u
    return u    

def fy(t,x,dx,u,px,pu,py): 
    y = np.array([x[0]+x[1],x[0]-x[1]], np.float64)
    #y = np.array([x[0]], np.float64)
    #print "u is ",u,x,y
    return y   
    
t1 = 30.0
dt = 0.05
t = np.arange(t0,t1,dt)
nt = len(t)
pu = np.array([0.1], np.float64)
py = np.array([1], np.float64)

nx = 2
nu = 1
ny = 2#1

npx = 3
npu = 0
npy = 0

xv_Euler,yv_Euler = ode.eulerintegration(ode.fdxsample2,fu,fy,x0,t ,px,pu,py,nu,ny) 

xv_RK,yv_RK       = ode.rk4(ode.fdxsample2,fu,fy,x0,t,px,pu,py,nu,ny) 

if plotsOn:
    print "plotting results"
    #dat = np.loadtxt('TY.dat')
    #pl.plot(dat[:,1],dat[:,2],'r.-')

    pl.figure(figsize=(8,8))
    pl.plot(xv_Euler[:,0],xv_Euler[:,1],'b')
    pl.plot(xv_RK[:,0],xv_RK[:,1],'g')
    pl.legend(("F2PY-python-Euler","F2PY-python-RungeKutta-4"))
    pl.savefig('phaseSpace_Euler_RK.pdf')
    
    
    pl.figure(figsize=(8,8))
    pl.plot(t,xv_Euler[:,0],'b',t,xv_Euler[:,1],'g',
            t,xv_RK[:,0],'b:',t,xv_RK[:,1],'g:')
    pl.legend(('$x_1$ Euler','$x_2$ Euler','$x_1$ RK','$x_2$ RK'))    
    pl.savefig('timeHistory_Euler_RK.pdf')
    
if plotsOn:
    pl.figure(figsize=(8,8))
    pl.plot(t, xv_Euler[:,0], t, xv_RK[:,0])
    pl.savefig('timeHistory1stState_Euler_RK.pdf')
    
    pl.figure(figsize=(8,8))
    pl.plot(t, xv_Euler[:,1], t, xv_RK[:,1])
    pl.savefig('timeHistory2ndtState_Euler_RK.pdf')