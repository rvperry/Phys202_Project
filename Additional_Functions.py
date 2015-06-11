#%matplotlib inline
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
from IPython.html.widgets import interact, fixed

gamma=4.5e-8
py = np.linspace(-100,100,100)
r0 = 25.0
px = -py**2/(4.0*r0) + r0

def star_ic2(M,radii,n):
    ic=[]
    for i in range(len(radii)):
        r=radii[i]
        N=n[i]
        theta=np.linspace(0,2*np.pi,N)
        for j in range(N):
            x=r*np.cos(theta[j])
            y=r*np.sin(theta[j])
            v=np.sqrt(gamma*M/r)
            dx=v*np.sin(theta[j])
            dy=-v*np.cos(theta[j])
            ic.append(np.array([x,y,dx,dy]))
    return ic

def S_ic2(y,r0,M,S,verr,aerr):
    x=-y**2/(4*r0)+r0
    v=verr*np.sqrt(2.0*gamma*(M+S)/np.sqrt(x**2+y**2))
    a=aerr*np.arctan(2*r0/y)
    vx=v*np.cos(a)
    vy=-v*np.sin(a)
    return np.array([x,y,vx,vy])

def Star_ic(M,radii,n,b):
    ic=[]
    for i in range(len(radii)):
        r=radii[i]
        N=n[i]
        theta=np.linspace(0,2*np.pi,N)
        for j in range(N):
            x=r*np.cos(theta[j])+b[0]
            y=r*np.sin(theta[j])+b[1]
            v=np.sqrt(gamma*M/r)
            dx=v*np.sin(theta[j])
            dy=-v*np.cos(theta[j])
            ic.append(np.array([x,y,dx,dy]))
    return ic