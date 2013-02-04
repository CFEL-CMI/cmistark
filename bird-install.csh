#!/usr/bin/env tcsh
# install jkext in current user's home directory
# To use this installation
setenv PYTHONPATH /afs/desy.de/group/cfel/cfeld-cmi/bird/lib/python
python setup.py install --home=/afs/desy.de/group/cfel/cfeld-cmi/bird
