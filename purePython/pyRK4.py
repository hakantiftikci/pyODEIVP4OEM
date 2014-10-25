# -*- coding: utf-8 -*-
"""
Created on Mon Mar 17 17:20:35 2014

@author: tai
"""

def stateDerivative(t,x,dx,fu,fdx,px,pu):   
    fu(t,x,u,pu)
    fdx(t,x,u,px,dx)


def RK4(fdx,fu,fy,x0,tv,xv,yv,px,pu,py,nt):

    a = [ [0.0, 0.0, 0.0, 0.0] ,  [1.0/2.0,0.0,0.0,0.0] , [0.0,1.0/2.0,0.0,0.0] , [0.0,0.0,1.0,0.0] ]
    b = [ 1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0 ]
    c = [ 0.0, 1.0/2.0, 1.0/2.0, 1.0 ]


    t = tv(1)
    x = x0

    fu(t,x,u,pu)
    fdx(t,x,u,px,dx)
    fy(t,x,dx,u,y,px,pu,py)

    xv[1,:] = x
    yv[1,:] = y

    for i in range(nt):
        dt = tv(i+1)-tv(i)

        deltax = 0
        for i in range(4):
            xstage = x
            for j in range(i-1):
                xstage = xstage+dt*a[i,j]*k[:,j]
            stateDerivative(t+c[i]*dt,xstage,k[:,i],fu,fdx,px,pu)
            deltax = deltax + dt*b[i]*k[:,i]
   
        x = x + deltax
        t = t + dt
        
        fu(t,x,u,pu)
        fdx(t,x,u,px,dx)
        fy(t,x,dx,u,y,px,pu,py)

        xv[1+m,:] = x
        yv[1+m,:] = y

