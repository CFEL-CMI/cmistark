#!/usr/bin/env tcsh
# install jkext in current user's home directory
# To use this installation, run the same setenv at login or before using CMIstark
setenv PYTHONUSERBASE $HOME/.python
python setup.py install --user
