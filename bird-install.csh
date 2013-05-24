#!/usr/bin/env tcsh
# install jkext for CFEL-CMI use on BIRD (DESY compte cluster)
setenv PYTHONPATH /afs/desy.de/group/cfel/cfeld-cmi/bird/lib/python
python setup.py install --home=/afs/desy.de/group/cfel/cfeld-cmi/bird
