#!/usr/bin/env python
# -*- coding: utf-8; fill-column: 120 -*-

from numpy import cos, exp, mgrid, pi, sin, zeros
from scipy.special import sph_harm

from enthought.mayavi import mlab
import jkext.wigner
from cmath import phase
# Create a sphere
r = 0.3
phi, theta = mgrid[0:2*pi:101j, 0:pi:101j]

x = r*sin(theta)*cos(phi)
y = r*sin(theta)*sin(phi)
z = r*cos(theta)
J=2 
K=1

for M in range(J+1):
    for K in range(0,2):
        if K == 0:
            srange = [0]
        else:
            srange = [0,1]
        for S in srange:
            fig = mlab.figure("asym J=%(J)d K=%(K)d M=%(M)d S=%(S)d" % {"J": J,"K": K, "M": M, "S": S}, bgcolor=(1, 1, 1), fgcolor=(0, 0, 0), size=(400, 300))
            fig.scene.disable_render = True
            mlab.clf()
            psi = zeros((101,101))
            psiphase = zeros((101,101))
            j = 0
            i = 0
            for dtheta in theta[0]:
                for dphi in phi[:,0]:
                    temp = jkext.wigner.drot(J, -K, -M, dtheta)*exp(-1j*K*dphi)+(-1)**(S)*jkext.wigner.drot(J, K,-M, dtheta)*exp(1j*K*dphi)
                    psi[i,j] = (temp.real**2 + temp.imag**2)
                    psiphase[i,j] = phase(temp)
                    i = i+1
                j = j+1
                i = 0
            mlab.mesh(x*psi,y*psi,z*psi,scalars=psiphase,vmin=-pi,vmax=pi)
            # phase will be between -pi and pi from blue to red. I.e we will have a discontinuity at changing from pi to - pi in color
            fig.scene.disable_render = False
            



mlab.show()
