Installing CMIstark
======================

Prerequisites and obtaining CMIstark
---------------------------------------

Since CMIstark is written in Python, you need to install Python; CMIstark requires Python version
3.4 or higher.

In addition, you need various Python extension packages, these include

* NumPy
* SciPy
* Tables
* Matplotlib

The CMIstark source code is currently only avaiable at https://stash.desy.de, please contact Jochen
Küpper <jochen.kuepper@cfel.de> for furhter details.


Installing CMIstark
----------------------

A normal installation is performed by simply running the command::

  python setup.py install

However, often you do not have the administrative rights to install in global directories, or simply
do not want to overrride a global installtion. In this case, you might want to perform a local
installation in your user directory using::

  python setup.py install --user

A similar setup can be achieved using::

  python setup.py develop --user

which, however, sets up the installation in such a way that changes to your source directory are
automatically and immediately visible through the installed version. This avoids repeated
re-installs while you are developing code.

Once you are satisfied with your changes you might consider reinstalling using one of the above two
options.

Fur further details of ``develop`` install, see https://pythonhosted.org/setuptools/setuptools.html#development-mode

Installing CMIstark: in user-specified path 
-------------------------------------------

Use PYTHONUSERBASE to specify the installation path::

  setenv PYTHONUSERBASE $HOME/.local
  python setup.py install --user

In the above example of installation (in tcsh shell), the module will be installed in the following path::

  $HOME/.local/lib/python/site-packages

and the scripts will be installed in the following path::

  $HOME/.local/bin

To import modules and call scripts of such user-specific installation, the following environment 
declarifications are required::

  setenv PATH /opt/local/bin:$HOME/.local/bin:$PATH
  setenv PYTHONUSERBASE $HOME/.local

The above example is provided for the tcsh shell. You can also then use ``site`` module of python
in python command prompt to make sure the environment is properly set up. For example::

  >>> import site
  >>> site.USER_BASE
  '$HOME/.local'

Also type "which ``name of script file``" to find the real path of the script called. It should
be in "$HOME/.local/bin".

For further details, see https://docs.python.org/3/install/index.html#inst-alt-install-user and
https://docs.python.org/3/using/cmdline.html#envvar-PYTHONUSERBASE  

.. comment
   Local Variables:
   coding: utf-8
   fill-column: 100
   truncate-lines: t
   End:
