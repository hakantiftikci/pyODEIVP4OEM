# -*- coding: utf-8 -*-
"""
Created on Mon Nov 24 14:42:56 2014

@author: tai
"""

import numpy as np
import pylab as pl   # use donly in test part 
import pyODEIVPSPV as ode

class Experiment:
    def __init__(self,fu,t,z,x0,p0,system):
        self.t = t
        self.z = z
        self.fu = fu
        self.x0 = x0
        self.pest = p0
        self.system = system
        
    def covarianceEstimate(self,p):                
        sim = self.system.generateVirtualExperiment(self.x0, self.fu, self.t, p)

        nz = self.z.shape[1]        
        nt = self.z.shape[0]
        R = np.zeros((nz,nz), np.float32)  # initialize covariance matrix to 0 for incremental computation
        for i in range(nt): 
            # for all measurements do
            delta = self.z[i,:] - sim.z[i,:]
            R = R + np.outer(delta,delta)
        
        return R        
        
    def computeFG(self,p,R,abseps,releps):        
        
        RI = np.linalg.inv(R)
        
        simnom = self.system.generateVirtualExperiment(self.x0, self.fu, self.t, p)
        
        _np = self.system.np
        _nt = self.z.shape[0]
        _ny = self.z.shape[1]
        
        
        Y = np.zeros( (_nt, _ny, len(p)), np.float32)
        
        for i in range(len(p)):        
            dp = 0*p
            
            if (p[i]==0):
                h = abseps                
            else:
                h = releps*p[i]
                
            dp[i] = h
                
            ppert = p.copy() # copy current thetax (before perturbation)
            ppert += dp
                
            simpert = self.system.generateVirtualExperiment(self.x0, self.fu, self.t, ppert)
            
            Y[:,:,i] = (simpert.z-simnom.z)/h

        F = np.zeros( (_np,_np), np.float32)            
        G = np.zeros( (_np,), np.float32)            
        for i in range(_nt):
            delta = self.z[i,:] - simnom.z[i,:]
            ytheta = Y[i,:,:]
            
            #print "ytheta size", ytheta.shape
            #print "delta  size", delta.shape
            #print "RI     size", RI.shape
            
            F += np.dot(ytheta.T , np.dot(RI,ytheta))
            G += np.dot(ytheta.T , np.dot(RI,delta))
            
            
        return (F,G,simnom,Y)
        
        
        
class System:
    def __init__(self,fdx,fy,p,nx,nu,ny,np,descr,fdx_extarg=(),fu_extarg=(),fy_extarg=()):
        self.fdx = fdx
        self.fy = fy
        self.p = p
        
        self.fdx_extarg = fdx_extarg
        self.fu_extarg  = fu_extarg
        self.fy_extarg  = fy_extarg

        # some sizes can be initialized from input vectors
        self.nx = nx
        self.nu = nu
        self.ny = ny       
        self.np = np

        #self.none = np.array([999.0], np.float32) # dummy array to pass as null array argument
        
        self.descr = descr
        
    def report(self):
        return ""
        
    def __str__(self):
        return "system:\n"+self.descr
        
        
        
    def generateVirtualExperiment(self, x0, fu, t, p, mu=(), sigma=()):        
        #if p==None:
        #p = self.p
                
        #xv,yv = rk4(fdx,fu,fy,x0,tv,p,nu,ny,[nx,np,nt,fdx_extra_args,fu_extra_args,fy_extra_args])
        nt = t.shape[0]
        xv,yv = ode.rk4(self.fdx, fu, self.fy, x0, t, p, self.nu, self.ny, self.nx, self.np, nt, self.fdx_extarg, self.fu_extarg, self.fy_extarg)
        
        z = yv.copy()
        for i in range(min(len(mu),len(sigma),yv.shape[1])):
            print "adding noise to %d output" % (i,)
            noise = np.random.normal(mu[i],sigma[i],z.shape[0])
            #print noise 
            z[:,i] += noise
            
        return Experiment(fu, t, z, x0, p, self)
        
    

class OutputError:
    def __init__(self,system, pest):
        self.system = system
        self.pest = pest
        self.experiments = []        
        
        
    def addExperiment(self,exp):
        self.experiments.append(exp)       
            
    def OEiters(self,niter):
        """ computes F,G matrices """
        p = self.pest.copy()
        cost = []
        parameters = []
        outputs = []
        sigmas = []
        for i in range(niter):
            exp = self.experiments[0]
            print "iteration ",i
            R = exp.covarianceEstimate(p)
            J = np.linalg.det(R)
            cost.append(J)
            parameters.append(p.copy())
            
            #,1e-6,1e-5
            #print R,np.linalg.det(R)
            
            F,G,znom,Y = exp.computeFG(p,R,1e-6,1e-5)
            invF = np.linalg.inv(F)
            deltheta = np.dot( invF, G)
            
            outputs.append(znom)
            sigmas.append(np.sqrt(np.diag(invF)))
            
            p  += deltheta
            #self.thetax_est  += ((G/np.sqrt(np.dot(G,G)))*0.1)
        return cost,parameters,outputs,sigmas


if __name__=="__main__":
    def fu1(t,x,pu):
        u = [10*np.sin(t)]
        #print "u is ",u
        return u    
    
    def fy(t,x,dx,u,px,pu,py): 
        y = np.array([x[0]+x[1],x[0]-x[1]], np.float32)
        return y   

    plotsOn = True
        
    nop = 3
    nox = 2
    nou = 1
    noy = 2
    p = np.array([1, 0.2, 3], np.float32)
    x01 = np.array([1,2], np.float32)
    pest = np.array([1.4, 0.4, 2.1], np.float32)
    pest = np.array([2.4, 0.4, 6.1], np.float32)
    # def system.__init__(self,fdx,fy,p,nx,nu,ny,np,descr):
    smd = System(ode.fdxsample2, fy, p, nox, nou, noy, nop, "SpringMassDamper system." )    
    
    tv1 = np.arange(0.0, 30.0, 0.1)
    exp1 = smd.generateVirtualExperiment(x01, fu1, tv1, p, (0.0, 0.0), (0.8,0.36))
    
    if plotsOn:
        pl.figure()
        pl.plot(exp1.z[:,0],exp1.z[:,1],'.-')
        
        pl.figure()
        pl.plot(tv1,exp1.z[:,0],'.-',tv1,exp1.z[:,1],'.-')    
        
    
    R = exp1.covarianceEstimate(pest)
    print "R = ", R
    
    abstol,reltol = 1e-6,1e-5
    F,G,simnom,Y = exp1.computeFG(pest,R,abstol,reltol)
    
    OE = OutputError(smd, pest)
    OE.addExperiment( exp1 )
    
    niter = 15
    costs,parameters,znominals,sigmas = OE.OEiters(niter)
    
    parameters = np.array(parameters)
    sigmas = np.array(sigmas)
    
    pl.figure()
    if False:
        pl.plot(range(niter),parameters[:,0],'+-',
                range(niter),parameters[:,1],'x-',
                range(niter),parameters[:,2],'^-')
        pl.errorbar(range(niter),parameters[:,0],sigmas[:,0])
        pl.errorbar(range(niter),parameters[:,1],sigmas[:,1])
        pl.errorbar(range(niter),parameters[:,2],sigmas[:,2])            
    else:
        if False:
            pl.plot(range(niter),parameters[:,0],'ro',mec='r')
            pl.plot(range(niter),parameters[:,1],'go',mec='g')
            pl.plot(range(niter),parameters[:,2],'bo',mec='b')
            pl.errorbar(range(niter),parameters[:,0],sigmas[:,0])
            pl.errorbar(range(niter),parameters[:,1],sigmas[:,1])
            pl.errorbar(range(niter),parameters[:,2],sigmas[:,2])                    
        else:
            pl.subplot(311);
            #pl.plot(range(niter),parameters[:,0],'ro',mec='r')
            pl.errorbar(range(niter),parameters[:,0],sigmas[:,0],fmt='r',ecolor='r')
            
            pl.subplot(312);
            #pl.plot(range(niter),parameters[:,1],'go',mec='g')
            pl.errorbar(range(niter),parameters[:,1],sigmas[:,1],fmt='g',ecolor='g')
            
            pl.subplot(313);
            #pl.plot(range(niter),parameters[:,2],'bo',mec='b')
            pl.errorbar(range(niter),parameters[:,2],sigmas[:,2],fmt='b',ecolor='b')                                
    
    pl.figure()
    pl.semilogy(range(niter),sigmas[:,0],'ro-',mec='r')
    pl.semilogy(range(niter),sigmas[:,1],'bo-',mec='b')
    pl.semilogy(range(niter),sigmas[:,2],'go-',mec='g')    
    
    pl.figure()
    pl.semilogy(range(niter),costs,'o-')
    pl.grid()    
    
    def colormap(i,niter): 
        l = float(i)/niter
        return  (l,1-l,0.0)
    

    pl.figure()
    for i in range(niter):
        pl.plot(znominals[i].t,znominals[i].z[:,0],color=colormap(i,niter));
    pl.plot(exp1.t,exp1.z[:,0],'.',mfc='b',mec='b',alpha=0.36,ms=8)        
    pl.grid()
    
    
    pl.figure()
    for i in range(niter):
        pl.plot(znominals[i].t,znominals[i].z[:,1],color=colormap(i,niter));
    pl.plot(exp1.t,exp1.z[:,1],'.',mfc='b',mec='b',alpha=0.36,ms=8)        
    pl.grid()

    pl.figure()
    for i in range(niter):
        pl.plot(znominals[i].z[:,0],znominals[i].z[:,1],color=colormap(i,niter));
    pl.plot(exp1.z[:,0],exp1.z[:,1],'.',mfc='b',mec='b',alpha=0.2,ms=3)            
    pl.grid()

