from __future__ import division

import numpy as num
import tables

def writeTable(file,tabletemplate,param,filters=tables.Filters(complevel=1, complib='zlib')):
    """
    ToDo doc
    """
    group = file.root
    try:
        table = file.getNode(group,'params')
    except tables.exceptions.NodeError:
        table = file.createTable(group, 'params', tabletemplate, "Calculation Parameters",filters=filters)
    data = table.row
    for x in tabletemplate.columns:
        data[x] = getattr(param,x)
        
    #data['type'] = param.type
    #data['Jmax_calc'] = param.Jmax_calc
    #data['Jmax_save'] = param.Jmax_save
    #data['isomer'] = param.isomer
    #data['saveevec'] = param.saveevec
    #data['rotcon'] = param.rotcon
    #data['quartic'] = param.quartic
    #data['dipole'] = param.dipole
    #data['polarizability'] = param.polarizability
    #data['watson'] = param.watson
    #data['symmetry'] = param.symmetry
    
    data.append()
    table.flush()


def readTable(file,tabletemplate,param):
    """
    ToDo doc
    """
    group = file.root
    name = "params"
    table = file.getNode(group,name)
    for x in tabletemplate.columns:
        setattr(param,x,table.col(x)[0]) # we just take the first element for now
        
    
