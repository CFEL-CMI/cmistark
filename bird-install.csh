#!/usr/bin/env tcsh
# install cmistark for CFEL-CMI use on BIRD (DESY compte cluster)
source /afs/desy.de/group/cfel/cfeld-cmi/bird/setup.csh
python setup.py install --home=${CMIBIRDPATH}
