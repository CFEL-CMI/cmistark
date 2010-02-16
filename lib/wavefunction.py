from __future__ import division

import numpy as num
import jkstark.wigner
import jkstark.starkeffect as starkeffect
from jkext.state import State

def asymwavefunction(J, K, M, S, theta):
    if K == 0:
        psi = (J+0.5)**(1/2)*jkstark.wigner.drot(J, -K, -M, theta) 
    else:
        psi = (J+0.5)**(1/2)*1/(2**(2/2))*(jkstark.wigner.drot(J, -K, -M, theta) +jkstark.wigner.drot(J, K, -M, theta))
    return psi



def calcwaveparam(state,dcfield,acfield,mol):
    
    param = starkeffect.CalculationParameter
    mol.getparam(param)
    acfield = self.testacfield(state,acfield)
    eigvectors = mol.geteigvectors(state,acfield,dcfield)
    Jmax = param.Jmax_save
    J = state.J()
    Ka = state.Ka()
    Kc = state.Kc()
    M = state.M()
    Jmin = M
    i = 0
    j = 0
    thetas = num.linspace(0,num.pi,201)
    psi = num.zeros_like(thetas)
    if param.symmetry == "C2a":
        if Ka%2 == 0:
            for theta in thetas:
                for J in range(Jmin,Jmax+1): # looping over the even states i.e E- and E+
                    Jeven = J-J%2
                    for K in range(Jeven,0,-2):
                        S = 1
                        psi[j] += asymwavefunction(J, K, M, S, theta)*eigvectors[i]
                        i = i +1
                    for K in range(0,Jeven+2,2):
                        S = 0
                        psi[j] += asymwavefunction(J, K, M, S, theta)*eigvectors[i]
                        i = i +1
                j = j+1
                i = 0
        else:
            psi = psi + 1
            #raise NotImplementedError("We only handle even Ka now")
    step = (1*(thetas < num.pi/2)+0.5*(thetas == num.pi/2))*psi
    integ = num.trapz(psi**2*num.sin(thetas),thetas)
    expcos = num.trapz(psi**2*num.cos(thetas)*num.sin(thetas),thetas)/integ
    exp2cos = num.trapz(psi**2*num.cos(thetas)**2*num.sin(thetas),thetas)/integ
    up_total = num.trapz(step**2*num.sin(thetas),thetas)/integ
    return expcos,exp2cos,up_total,psi,thetas
