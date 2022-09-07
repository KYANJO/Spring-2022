#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 25 19:02:33 2022

@author: Brian KYANJO
"""

'''
Parameters used:
----------------
ql  - Is an array containing height(hl) and momentum(hul) left states ([hl,hul])
qr  - Is an array containing height(hr) and momentum(hur) right states ([hr,hur])
qm  - Is an array containing height(hm) and momentum(hum) intermidiate states ([hm,hum])
qms - Is an array containing height(hms) and momentum(hums) intermidiate shock states ([hms,hums])
qmr - Is an array containing height(hmr) and momentum(humr) intermidiate rarefaction states ([hmr,humr])
g   - gravity
 
'''
g = 1

from numpy import *
import time

def pospart(x):
    '''
    Returns a value greater than zero
    '''
    return max(1e-15,x)

#All shock solution
def ff(hm,um,g,hl,hr,ul,ur):
    '''
    Input:
    ------
    hm - Intemediate hieght field
    um - Intemediate velocity field
    Returns:
    ------
    f1 - A curve that corresponds to a 2-shock (connects qm to qr)
    f2 - A curve that corresponds to a 1-shock (connects ql to qm)
    '''
    #hugoniot_locus
    f1 = (um - (ur + (hm-hr)*sqrt((g/2)*(1/hm + 1/pospart(hr))))) 
    f2 = (um - (ul - (hm-hl)*sqrt((g/2)*(1/hm + 1/pospart(hl)))))
    return f1,f2
    
#Derivative of function ff
def dff(hm,um,g,hl,hr,ul,ur):
    '''
    Input:
    ------
    hm - Intemediate hieght field
    um - Intemediate velocity field
    Returns:
    --------    
    f1h - Derivative of f1 wrt h
    f1u - Derivative of f1 wrt u
    f2h - Derivative of f2 wrt h
    f2u - Derivative of f2 wrt u
    '''
    f1h = (sqrt(2*g*(hm + hr)/(hm*pospart(hr))))*(-2*hm*(hm + hr) + hr*(hm - hr)) \
          / (4*hm*(hm + hr))
   
    f2h = (sqrt(2*g*(hm + hl)/(hm*pospart(hl))))*(2*hm*(hm + hl) + hl*(hl - hm)) \
          / (4*hm*(hm + hl))
          
    f1u = 1
    f2u = 1
    
    return f1h,f1u,f2h,f2u

#Jacobian
def J(dff,hm,um,g,hl,hr,ul,ur):
    '''
    Input:
    ------
    hm - Intemediate hieght field
    um - Intemediate velocity field
    Returns: A 2x2 matrix that contains derivatives of function ff wrt h and u.
    -------
    '''
    f1h,f1u,f2h,f2u = dff(hm,um,g,hl,hr,ul,ur)
    
    return array([[f1h,f1u],[f2h,f2u]])

#Inverse of Jacobian transpose and Jacobian
def Jinv(hm,um,hl,hr,ul,ur,g):
    '''
    Input:
    ------
    hm - Intemediate hieght field
    um - Intemediate velocity field
    Returns: The inverse of the Jacobian transpose and Jacobian matrix multiplied by the Jacobian
    -------
    '''
    JT = J(dff,hm,um,g,hl,hr,ul,ur).T@J(dff,hm,um,g,hl,hr,ul,ur) 
    return linalg.inv(JT)@J(dff,hm,um,g,hl,hr,ul,ur).T

def f(hm,um,g,hl,hr,ul,ur):
    '''
    Input:
    ------
    hm - Intemediate hieght field
    um - Intemediate velocity field
    Returns: An array of all shock solutions (f1 - 2-shock and f2 - 1-shock)
    -------
    '''
    f1,f2 = ff(hm,um,g,hl,hr,ul,ur)
    return array([f1,f2])

#shock wave solution
def GN(ql,qr,g):
    '''
    Description: Newton solver used to generate all shock Riemann solution
    -----------
    Returns:
    -------
    hm - Intermediate shock hieght field
    um - Intermediate shock velocity field
    '''
    #max _iterations
    max_iter = 20
    #tolerance
    epsilon  = 1e-16
    
    #intial conditions (IC)
    ho = 0.1
    uo = 0.01
    
    #left state intial height and momentum field 
    hl = ql[0]
    hul = ql[1]
    
    #right state intial height and momentum field
    hr = qr[0]
    hur = qr[1]
    
    #left and right intial states velocities
    ul = hul/pospart(hl)
    ur = hur/pospart(hr)

    #intial value
    vo = array([ho,uo])

    #Newton solver
    for i in range(max_iter):
        
        dq = -Jinv(ho,uo,hl,hr,ul,ur,g)@f(ho,uo,g,hl,hr,ul,ur)

        v1 = vo + dq

        if linalg.norm(v1-vo) < epsilon:
            break
        else:
            vo = v1
            ho = v1[0]
            uo = v1[1]
    
    #Intermediate fields
    hm = v1[0]
    um =v1[1]
    
    return hm,um

#Riemann Solvers:
def GN_solver(Q_ext,GN):
    """  Input : 
            Q_ext : Array of N+4 Q values.   Boundary conditions are included.
            
        Output : 
            waves  : Jump in Q at edges -3/2, -1/2, ..., N-1/2, N+1/2 (N+3 values total)
            speeds : Array of speeds (N+3 values)
            apdq   : Positive fluctuations (N+3 values)
            amdq   : Negative fluctuations (N+3 values)
        """    
        
     # jump in Q at each interface
    delta = Q_ext[1:,:]-Q_ext[:-1,:]

    d0 = delta[:,[0]]
    d1 = delta[:,[1]]
    
    qold1 = Q_ext[:,0]
    qold2 = Q_ext[:,1]
    
    mx = delta.shape[0]
    
    # Array of wave 1 and 2
    w1 = zeros(delta.shape)
    w2 = zeros(delta.shape)
    
    # Array of speed 1 and 2
    s1 = zeros((delta.shape[0],1))
    s2 = zeros((delta.shape[0],1))
    
    amdq = zeros(delta.shape)
    apdq = zeros(delta.shape)
    
    for i in range(1,mx):
        
        ql = array([qold1[i],qold2[i]])
        qr = array([qold1[i+1],qold2[i+1]]) #at edges
            
        #at the intefaces
        hms,ums = GN(ql,qr,g)
        hums = hms*ums
        
        #state at the interface
        qm = array([hms,hums])
        
        #fluctuations
        amdq[i] = flux(qm) - flux(ql)
        apdq[i] = flux(qr) - flux(qm)
        
    return amdq,apdq



def forestclaw_solver(Q_ext,meqn):
    """  Input : 
            Q_ext : Array of N+4 Q values.   Boundary conditions are included.
            
        Output : 
            waves  : Jump in Q at edges -3/2, -1/2, ..., N-1/2, N+1/2 (N+3 values total)
            speeds : Array of speeds (N+3 values)
            apdq   : Positive fluctuations (N+3 values)
            amdq   : Negative fluctuations (N+3 values)
    """    
        
    # jump in Q at each interface
    delta = Q_ext[1:,:]-Q_ext[:-1,:]
    
    # For most problems, the number of waves is equal to the number of equations
    mwaves = meqn

    d0 = delta[:,[0]]
    d1 = delta[:,[1]]
    
    h = Q_ext[:,0]
    u = Q_ext[:,1]/Q_ext[:,0]
    n = delta.shape[0]
    
    # Array of wave 1 and 2
    w1 = zeros(delta.shape)
    w2 = zeros(delta.shape)
  
    # Array of speed 1 and 2
    s1 = zeros((delta.shape[0],1))
    s2 = zeros((delta.shape[0],1))
   
    for i in range(1,n):
        u_hat = (sqrt(h[i-1])*u[i-1]+sqrt(h[i])*u[i])/(sqrt(h[i-1])+sqrt(h[i]))
        h_bar = (1/2)*(h[i-1]+h[i])
        c_hat = sqrt(g*h_bar) 
        
        # Eigenvalues
        l1 = u_hat - c_hat        
        l2 = u_hat + c_hat   

        # Eigenvectors
        r1 = array([1, l1])       
        r2 = array([1, l2])          
        
        R = array([r1,r2]).T
        
        # Vector of eigenvalues
        evals =  array([l1,l2])     

        # Solve R*alpha = delta to get a1=alpha[0], a2=alpha[1]
        a1 = ((u_hat+c_hat)*d0[i-1]-d1[i-1])/(2*c_hat)
        a2 = (-(u_hat-c_hat)*d0[i-1]+d1[i-1])/(2*c_hat)
        
        # Wave and speed 1
        w1[i-1] = a1*R[:,[0]].T
        s1[i-1] = evals[0]

        # Wave and speed 2
        w2[i-1] = a2*R[:,[1]].T
        s2[i-1] = evals[1]
    
    waves = (w1,w2)             # P^th wave at each interface
    speeds = (s1,s2)            # Speeds at each interface

    # Fluctuations
    amdq = zeros(delta.shape)
    apdq = zeros(delta.shape)
    for mw in range(mwaves):
        sm = where(speeds[mw] < 0, speeds[mw], 0)
        amdq += sm*waves[mw]
        
        sp = where(speeds[mw] > 0, speeds[mw], 0)
        apdq += sp*waves[mw]
    
    return amdq,apdq

#solver based on clawpack
def claw(ax, bx, mx, Tfinal, umax, ql, qr, cfl,h,
          meqn=1, \
          switch=1,\
          GN=None,\
          qinit=None, \
          solver=None, 
          bc=None):

    dx = (bx-ax)/mx
    xe = linspace(ax,bx,mx+1)  # Edge locations
    xc = xe[:-1] + dx/2       # Cell-center locations
    
    dt_est = cfl*dx/umax;
    nout = int(floor(Tfinal/dt_est) + 1)

    # For many problems we can assume mwaves=meqn.
    mwaves = meqn
        
    # Temporal mesh
    t0 = 0
    tvec = linspace(t0,Tfinal,nout+1)
    dt = Tfinal/nout
    
    assert qinit is not None, 'No user supplied initialization routine'
    assert bc is not None,    'No user supplied boundary conditions'

    # Initial the solution
    ej = array([1,1])
    if switch == 0:
        q0 = qinit(xc,meqn,ql+h*ej,qr+h*ej)    # Should be [size(xc), meqn]
    elif switch == 1:
        q0 = qinit(xc,meqn,ql-h*ej,qr-h*ej)    # Should be [size(xc), meqn]
    else:
        q0 = qinit(xc,meqn,ql,qr)    # Should be [size(xc), meqn]
        
    # Store time stolutions
    Q = empty((mx,meqn,nout+1))    # Doesn't include ghost cells
    Q[:,:,0] = q0

    q = q0
    
    dtdx = dt/dx
    t = t0
    for n in range(0,nout):
        t = tvec[n]
        
        # Add 2 ghost cells at each end of the domain;  
        q_ext = bc(q)

        # Get fluctuations
        if solver == 0:
            start_GN = time.time()
            amdq, apdq = GN_solver(q_ext,GN)
            end_GN = time.time()
        elif solver == 1:
            start_fs = time.time()
            amdq, apdq = forestclaw_solver(q_ext,meqn)
            end_fs = time.time()
            
        # First order update
        q = q - dtdx*(apdq[1:-2,:] + amdq[2:-1,:])
        
        Q[:,:,n+1] = q
    if solver == 0:
        #print("\n")
        print('solver used is: forestclaw')
        print("dt = {:.4e}".format(dt))
        print("Number of time steps = {}".format(nout))
        print("Computational time = {:.4e}".format(end_GN - start_GN),"seconds")
            
    elif solver == 1:
        print('solver used is: Gauss_Newton based')
        print("dt = {:.4e}".format(dt))
        print("Number of time steps = {}".format(nout))
        print("Computational time = {:.4e}".format(end_fs - start_fs),"seconds")
    return Q, xc, tvec


def h_init(x,hl,hr):    
    q0 = where(x < 0,hl,hr)
    return q0

def hu_init(x,hl,ul,hr,ur):    
    #q0 = zeros(x.shape)  
    q0 = where(x<0,hl*ul,hr*ur)
    return q0

def qinit(x,meqn,ql,qr):
    #initial height fields(left and right)
    hl = ql[0]
    hr = qr[0]
    
    #initial momentum fields(left and right)
    hul = ql[1]
    hur = qr[1]
    
    #initial momentum fields(left and right)
    ul = hul/pospart(hl)
    ur = hur/pospart(hr)

    q = zeros((x.shape[0],meqn))
    q[:,0] = h_init(x,hl,hr)
    q[:,1] = hu_init(x,hl,ul,hr,ur)
    
    return q

# Boundary conditions
def bc_extrap(Q):
    """ Extend Q with extrapolation boundary conditions """
        
    Q_ext = concatenate((Q[[1,0],:], Q, Q[[-1,-2],:]))
    return Q_ext

#flux
def flux(q):
    '''
    input:
    -----
    q - state at the interface
    return:
    -------
    f - flux at the interface
    '''
    q1 = q[0]
    q2 = q[1]
    f = zeros(2)
    f[0] = q2
    f[1] = (((q2)**2)/q1) + (0.5*g*(q1)**2)
    return f


def pospart(x):
    '''
    Returns a value greater than zero
    '''
    return max(1e-15,x)


