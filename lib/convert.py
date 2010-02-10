from __future__ import division

"""Various conversions for easy comparison with F and H articles"""
import numpy
import jkext.const

def dcfields2omega(dcfields, rotcon, dipole):
    """ToDo: document"""
    return numpy.array(dcfields)*dipole/rotcon

def acfields2deltaomega(acfields, rotcon, polarizability):
    """ToDo: document"""
    return acfields**2*polarizability/(4*rotcon)    

def omega2dcfields(omega, rotcon, dipole):
    """ToDo: document"""
    assert 0 != dipole
    return numpy.array(omega)*rotcon/dipole

def deltaomega2acfields(deltaomega, rotcon, polarizability):
    """ToDo: document"""
    assert 0 != polarizability
    return (numpy.array(deltaomega)*4*rotcon/polarizability)**0.5

def w_cm22V_m(val):
    """Convert intensity in w/cm2 to field"""
    return (1e4*val*2*jkext.const.vacuum_permeability*jkext.const.speed_of_light)**(1/2)

def V_m2w_cm2(val):
    """Convert intensity in w/cm2 to field"""
    return (val**2/(1e4*2*jkext.const.vacuum_permeability*jkext.const.speed_of_light))

