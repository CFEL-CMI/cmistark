from numpy import cos, exp, mgrid, pi, sin, zeros
from scipy.special import sph_harm

from enthought.mayavi import mlab

import jkext.wigner

# Create a sphere
r = 0.3
phi, theta = mgrid[0:2*pi:101j, 0:pi:101j]

x = r*sin(theta)*cos(phi)
y = r*sin(theta)*sin(phi)
z = r*cos(theta)

J=2
K=1

plot_gsl = False
if K == 0:
    plot_sh = True
else:
    plot_sh = False

for M in range(J+1):
    fig = mlab.figure(100*J+10*K+M, bgcolor=(1, 1, 1), fgcolor=(0, 0, 0), size=(400, 300))
    fig.scene.disable_render = True
    mlab.clf()
    s = zeros((101,101))
    j = 0
    i = 0
    for dtheta in theta[0]:
        s[:,j] = (jkext.wigner.drot(J, K, M, dtheta))**2
        j = j+1
    mlab.mesh(s*x,s*y,s*z)
    fig.scene.disable_render = False

    if plot_sh == True:
        fig = mlab.figure(1000+100*J+10*K+M, bgcolor=(1, 1, 1), fgcolor=(0, 0, 0), size=(400, 300))
        mlab.clf()
        fig.scene.disable_render = True
        s = sph_harm(M, J, phi, theta)
        t = (s.real**2 + s.imag**2) # abs(s)**2
        mlab.mesh(t*x, t*y, t*z)
        fig.scene.disable_render = False

    if plot_gsl == True:
        fig.scene.disable_render = True
        fig = mlab.figure(2, bgcolor=(1, 1, 1), fgcolor=(0, 0, 0), size=(400, 300))
        mlab.clf()
        s=zeros((101,101))
        j=0
        i=0
        for dtheta in theta[0]:
            i=0
            s[:,j]=jkext.wigner.drot_gsl(J,K,M,dtheta)
            j=j+1
        mlab.mesh(s*x,s*y,s*z)
        fig.scene.disable_render = False

mlab.show()
