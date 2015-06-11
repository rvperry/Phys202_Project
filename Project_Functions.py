#%matplotlib inline
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
from IPython.html.widgets import interact, fixed
import Project_Functions as pf

gamma=4.5e-8 #kpc^3/(SolarMass)(age^2) age=10^8 years
py = np.linspace(-100,100,100)
r0 = 25.0
px = -py**2/(4.0*r0) + r0

def derivs(Rrv,t,M,S):
    rx,ry,vx,vy,Rx,Ry,Vx,Vy=tuple(Rrv)

    r=np.sqrt(rx**2+ry**2)
    R=np.sqrt(Rx**2+Ry**2)
    rhox=Rx-rx
    rhoy=Ry-ry
    rho=np.sqrt(rhox**2+rhoy**2)
 
    drx = vx
    dry = vy
    dvx = -gamma*(M*rx/(r**3)-S*rhox/(rho**3)+S*Rx/(R**3))
    dvy = -gamma*(M*ry/(r**3)-S*rhoy/(rho**3)+S*Ry/(R**3))
    dRx = Vx
    dRy = Vy
    dVx = -gamma*(M+S)*Rx/(R**3)
    dVy = -gamma*(M+S)*Ry/(R**3)
    return np.array([drx,dry,dvx,dvy,dRx,dRy,dVx,dVy])
    
def star_ic(M,radii,n):
    ic=[]
    for i in range(len(radii)):
        r=radii[i]
        N=n[i]
        theta=np.linspace(0,2*np.pi,N)
        for j in range(N):
            x=r*np.cos(theta[j])
            y=r*np.sin(theta[j])
            v=np.sqrt(gamma*M/r)
            dx=-v*np.sin(theta[j])
            dy=v*np.cos(theta[j])
            ic.append(np.array([x,y,dx,dy]))
    return ic
    
def S_ic(y,r0,M,S):
    x=-y**2/(4*r0)+r0
    v=np.sqrt(2.0*gamma*(M+S)/np.sqrt(x**2+y**2))
    a=np.arctan(2*r0/y)
    vx=v*np.cos(a)
    vy=-v*np.sin(a)
    return np.array([x,y,vx,vy])
    
def solve_one_star(ics, icg, M, S, tmax, ntimes):
    t = np.linspace(0, tmax, ntimes)
    ic = np.hstack([ics, icg])
    soln = odeint(derivs, ic, t, args=(M,S))
    return soln
    

    
def solve_all_stars(ics, icg, M, S, tmax, ntimes):
    solns = []
    for ic in ics:
        soln = solve_one_star(ic, icg, M, S, tmax, ntimes)
        solns.append(soln)
    return solns
    
def plot_one_soln(soln, j, lim):
    py = np.linspace(-100,100,100)
    r0 = 25.0
    px = -py**2/(4.0*r0) + r0
    plt.scatter(soln[j,0], soln[j,1])
    plt.scatter(soln[j,4], soln[j,5])
    plt.scatter(0,0)
    plt.plot(px,py)
    plt.xlim(-lim,lim)
    plt.ylim(-lim,lim)
    
def plot_all_solns(solns, j, lim):
    py = np.linspace(-100,100,100)
    r0=25.0
    px = -py**2/(4.0*r0) + r0
    x = np.array([soln[j,0] for soln in solns])
    y = np.array([soln[j,1] for soln in solns])
    X = np.array([soln[j,4] for soln in solns])
    Y = np.array([soln[j,5] for soln in solns])
    plt.scatter(x, y,color='green')
    plt.scatter(X, Y,color='gold')
    plt.scatter(0,0)
    plt.plot(px,py)
    plt.xlim(-lim,lim)
    plt.ylim(-lim,lim);
    
def subplot(solns,j,lim):    
    plt.figure(figsize=(15,7))
    plt.subplot(2,4,1)
    pf.plot_all_solns(solns,j[0],lim)
    plt.subplot(2,4,2)
    pf.plot_all_solns(solns,j[1],lim)
    plt.subplot(2,4,3)
    pf.plot_all_solns(solns,j[2],lim)
    plt.subplot(2,4,4)
    pf.plot_all_solns(solns,j[3],lim)
    plt.subplot(2,4,5)
    pf.plot_all_solns(solns,j[4],lim)
    plt.subplot(2,4,6)
    pf.plot_all_solns(solns,j[5],lim)
    plt.subplot(2,4,7)
    pf.plot_all_solns(solns,j[6],lim)
    plt.subplot(2,4,8)
    pf.plot_all_solns(solns,j[7],lim)