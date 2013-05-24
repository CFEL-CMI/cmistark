#!/usr/bin/env tcsh
# install jkext in current user's home directory
# To use this installation, run the same setenv at login or before using CMIstark
setenv PYTHONPATH $HOME/.python/lib/python
python setup.py install --home=$HOME/.python
