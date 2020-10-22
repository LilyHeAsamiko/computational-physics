# -*- coding: utf-8 -*-
"""
Created on Fri Feb 22 05:56:32 2019

@author: user
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from numpy import linalg as LA
from numpy import *


"""
Add basis functions p1,p2,q1,q2 here
"""


def init_1d_spline(x,f,h):
    # now using complete boundary conditions
    # with forward/backward derivative
    # - natural boundary conditions commented
    a=zeros((len(x),))
    b=zeros((len(x),))
    c=zeros((len(x),))
    d=zeros((len(x),))
    fx=zeros((len(x),))

    # a[0]=1.0 # not needed
    b[0]=1.0

    # natural boundary conditions 
    #c[0]=0.5
    #d[0]=1.5*(f[1]-f[0])/(x[1]-x[0])

    # complete boundary conditions
    c[0]=0.0
    d[0]=(f[1]-f[0])/(x[1]-x[0])
    
    for i in range(1,len(x)-1):
        d[i]=6.0*(h[i]/h[i-1]-h[i-1]/h[i])*f[i]-6.0*h[i]/h[i-1]*f[i-1]+6.0*h[i-1]/h[i]*f[i+1]
        a[i]=2.0*h[i]
        b[i]=4.0*(h[i]+h[i-1])
        c[i]=2.0*h[i-1]        
    #end for

    
    b[-1]=1.0
    #c[-1]=1.0 # not needed

    # natural boundary conditions
    #a[-1]=0.5
    #d[-1]=1.5*(f[-1]-f[-2])/(x[-1]-x[-2])

    # complete boundary conditions
    a[-1]=0.0
    d[-1]=(f[-1]-f[-2])/(x[-1]-x[-2])
    
    # solve tridiagonal eq. A*f=d
    c[0]=c[0]/b[0]
    d[0]=d[0]/b[0]
    for i in range(1,len(x)-1):
        temp=b[i]-c[i-1]*a[i]
        c[i]=c[i]/temp
        d[i]=(d[i]-d[i-1]*a[i])/temp
    #end for
        
    fx[-1]=d[-1]
    for i in range(len(x)-2,-1,-1):
        fx[i]=d[i]-c[i]*fx[i+1]
    #end for
        
    return fx
# end function init_1d_spline

""" 
Add smoothing functions 

def smooth1d(x,f,factor=3):
    ...
    ...
    return ...

def smooth2d(x,y,f,factor=3):
    ...
    ...
    return ... 

def smooth3d(x,y,z,f,factor=3):
    ...
    ...
    ...
    return ...
"""

class spline(object):

    def __init__(self,*args,**kwargs):
        self.dims=kwargs['dims']
        if (self.dims==1):
            self.x=kwargs['x']
            self.f=kwargs['f']
            self.hx=np.diff(self.x)
            self.fx=init_1d_spline(self.x,self.f,self.hx)
        elif (self.dims==2):
            self.x=kwargs['x']
            self.y=kwargs['y']
            self.f=kwargs['f']
            self.hx=np.diff(self.x)
            self.hy=np.diff(self.y)
            self.fx=zeros(shape(self.f))
            self.fy=zeros(shape(self.f))
            self.fxy=zeros(shape(self.f))
            for i in range(max([len(self.x),len(self.y)])):
                if (i<len(self.y)):
                    self.fx[:,i]=init_1d_spline(self.x,self.f[:,i],self.hx)
                if (i<len(self.x)):
                    self.fy[i,:]=init_1d_spline(self.y,self.f[i,:],self.hy)
            #end for
            for i in range(len(self.y)):
                self.fxy[:,i]=init_1d_spline(self.x,self.fy[:,i],self.hx)
            #end for
        elif (self.dims==3):
            self.x=kwargs['x']
            self.y=kwargs['y']
            self.z=kwargs['z']
            self.f=kwargs['f']
            self.hx=np.diff(self.x)
            self.hy=np.diff(self.y)
            self.hz=np.diff(self.z)
            self.fx=zeros(shape(self.f))
            self.fy=zeros(shape(self.f))
            self.fz=zeros(shape(self.f))
            self.fxy=zeros(shape(self.f))
            self.fxz=zeros(shape(self.f))
            self.fyz=zeros(shape(self.f))
            self.fxyz=zeros(shape(self.f))
            for i in range(max([len(self.x),len(self.y),len(self.z)])):
                for j in range(max([len(self.x),len(self.y),len(self.z)])):
                    if (i<len(self.y) and j<len(self.z)):
                        self.fx[:,i,j]=init_1d_spline(self.x,self.f[:,i,j],self.hx)
                    if (i<len(self.x) and j<len(self.z)):
                        self.fy[i,:,j]=init_1d_spline(self.y,self.f[i,:,j],self.hy)
                    if (i<len(self.x) and j<len(self.y)):
                        self.fz[i,j,:]=init_1d_spline(self.z,self.f[i,j,:],self.hz)
            #end for
            for i in range(max([len(self.x),len(self.y),len(self.z)])):
                for j in range(max([len(self.x),len(self.y),len(self.z)])):
                    if (i<len(self.y) and j<len(self.z)):
                        self.fxy[:,i,j]=init_1d_spline(self.x,self.fy[:,i,j],self.hx)
                    if (i<len(self.y) and j<len(self.z)):
                        self.fxz[:,i,j]=init_1d_spline(self.x,self.fz[:,i,j],self.hx)
                    if (i<len(self.x) and j<len(self.z)):
                        self.fyz[i,:,j]=init_1d_spline(self.y,self.fz[i,:,j],self.hy)
            #end for
            for i in range(len(self.y)):
                for j in range(len(self.z)):
                    self.fxyz[:,i,j]=init_1d_spline(self.x,self.fyz[:,i,j],self.hx)
            #end for
        else:
            print('Either dims is missing or specific dims is not available')
        #end if

    def p1(self, t):
        return (1+2*t)*((t-1)**2)  
    
    def p2(self, t):
        return t**2*(3-2*t)
    def q1(self, t):
        return t*(t-1)**2   
    def q2(self, t):
        return t**2*(t-1)
            
    def eval1d(self,x):
        if isscalar(x):
            # make sure that x is array
            x=array([x])
        N=len(self.x)-1
        f=zeros((len(x),))
        ii=0
        for val in x:
            # round and find the closest integer for i
            i=floor(where(self.x<=val)[0][-1]).astype(int)
            # when the axis reaches the maximum of the x, set the value of the
            # f as the last value for interpolation
            if i==N:
                f[ii]=self.f[i]
            # calculated according to the Hermite cubic splines
            else:
                t=(val-self.x[i])/self.hx[i]
                f[ii]=self.f[i]*self.p1(t)+self.f[i+1]*self.p2(t)+self.hx[i]*(self.fx[i]*self.q1(t)+self.fx[i+1]*self.q2(t))
            ii+=1

        return f
    #end eval1d
    
    def test1d(self, x, y, tol):
        plt.figure()
        plt.plot(y)
        plt.title('error norm1')
        plt.savefig('spline 1D error.pdf',dpi=200)
        err = LA.norm(y,2)
        print('the norm 2 of the difference between interpolation and analytic function:')
        print(np.array(err))
        if err < tol:
            print('1D interpolate is OK')
        return err


#end class spline


    
def main():

    fig1d = plt.figure()
    ax1d = fig1d.add_subplot(111)
    plt.title('interpolation of f')
    # 1d example
    x=np.linspace(-1.5,1.5,100)
    # test with sin
    def fun(x):
        return x
    y=fun(x)
    spl1d=spline(x=x,f=y,dims=1)
    xx=np.linspace(-1.5,0.1,100)
    ax1d.plot(xx,spl1d.eval1d(xx))   
    ax1d.plot(x,y,'o',xx,xx,'r--')
    print('the interpolation: %s'% repr(spl1d.eval1d(xx)[-1]))
    plt.savefig('spline 1D.pdf',dpi=200)
    tol = 0.5
    err1 = spl1d.test1d(xx, spl1d.eval1d(xx)-xx, tol)    

    

    plt.show()
#end main
    
if __name__=="__main__":
    main()
