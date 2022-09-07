#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 27 20:47:45 2022

@author: brian kyanjo
"""

import numpy as np
import statistics
from Gauss_Newton import *
from matplotlib.pylab import *


def Rieman_solvers(hl,ul,hr,ur, mq, umax, g,ax,bx, ay,by,mx,h, \
                   to, Tfinal,meqn,cfl):
    #intial data
    ql = array([hl,hl*ul])
    qr = array([hr,hr*ur])
    mo = array([ql,qr])
    print(mo.T)


    #Based on Gauss-Newton
    QN,x,tvec = claw(ax, bx, mx, Tfinal, umax, ql, qr, cfl,h,
              meqn=2, \
              switch=3,\
              GN=GN,\
              qinit=qinit,
              solver = 0,
              bc=bc_extrap)

    #Based on forestclaw: Roe-solver
    Qf,x,tvec = claw(ax, bx, mx, Tfinal, umax, ql, qr, cfl,h,
              meqn=2, \
              switch=3,\
              GN=GN,\
              qinit=qinit,
              solver = 1,
              bc=bc_extrap)
        
    Qfma,x,tvec = claw(ax, bx, mx, Tfinal, umax, ql, qr, cfl,h,
              meqn=2, \
              switch=0,\
              GN=GN,\
              qinit=qinit,
              solver = 1,
              bc=bc_extrap)
        
    Qfmin,x,tvec = claw(ax, bx, mx, Tfinal, umax, ql, qr, cfl,h,
              meqn=2, \
              switch=1,\
              GN=GN,\
              qinit=qinit,
              solver = 1,
              bc=bc_extrap)


    fig = figure(1)
    clf()

    if mq == 0:
        tstr = 'Height : t = {:.4f}'
    else:
        tstr = 'Momentum : t = {:.4f}'

    # Convert list of lists to list
    def flat(t):
        return [i for j in t for i in j]

    htitle = title(tstr.format(0),fontsize=18)

    #initialise the approximate soln with Roe-solver
    q0 = QN[:,mq,0]
    hdr, = plot(x,q0,'r.',markersize=5,label='$GN_{solver}$')
    
    alpha = 1e-16
    n  = 2
    std = 1e-5 #fixed standard deviation
    #std = statistics.stdev(Qf[:,mq,0])
    ep = std*array(flat(np.random.randn(len(Qf),1))) #+ statistics.mean(Qf[:,mq,0])
    qf = Qf[:,mq,0] + ep
    hdf, = plot(x,qf,'b.',markersize=5,label='$fclaw_{solver}$')

    for i,t in enumerate(tvec):

        #Using Gauss Newton
        q = QN[:,mq,i]
        hdr.set_ydata(q)

        #based on Forestclaw: Roe solver (simulated data)
        #std = statistics.stdev(Qf[:,mq,i])
        ep = std*array(flat(np.random.randn(len(Qf),1))) #+ statistics.mean(Qf[:,mq,i])
        qf = Qf[:,mq,i] + ep
        hdf.set_ydata(qf)
        
        # Noise
        epd = std*array((np.random.randn(len(Qf[:,:,i]),2)))
        
        # simulated data
        d = Qf[:,:,i] + epd
        print(d)
        # Approximate Jacobian
        J = (Qfma[:,:,i]  - Qfmin[:,:,i])/(2*h)
        print(J.shape)
        #parameter estimates
        m = linalg.pinv(J.T@J + (alpha**2)*eye(n))@(J.T@(Qf[:,:,i] - d) + (alpha**2)*eye(n)*mo.T)
        print(m)

        xlabel('x',fontsize=16)
        ylabel('q(x,t)',fontsize=16)
        htitle.set_text(tstr.format(t))

        legend()

        ylim([ay,by])

        pause(0.1)
        savefig('/Users/mathadmin/Documents/phd-coursework')

        fig.canvas.draw()                  