from __future__ import division

"""Various conversions for easy comparison with F and H articles"""
import numpy

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
