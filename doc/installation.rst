Installing CMIdiffract
======================

Prerequisites and obtaining CMIdiffract
---------------------------------------

Since CMIstark is written in Python, you need to install Python; CMIstark requires Python version
3.4 or higher.

In addition, you need various Python extension packages, this includes
* NumPy
* Tables

The CMIstark source code is currently only avaiable at https://stash.desy.de, please contact Jochen
Küpper <jochen.kuepper@cfel.de> for furhter details.


Installing CMIdiffract
----------------------

A normal installation is performed by simply running the command::

  python setup.py install

However, often you do not have the administrative rights to install in global directories, or simply
do not want to overrride a global installtion. In this case, you might want to perform a local
installation in your user directory using::

  python setup.py install --user

A similar setup can be achieved using::

  python setup.py develop --user

which, however, setups up the installation in such a way that changes to your source directory are
automaticall and immediately visible through the installed version. This avoids repeated re-installs
while you are developing code.

Once you are satisfied with your changes you might consider to reinstall using one of first two
options.



.. comment
   Local Variables:
   coding: utf-8
   fill-column: 100
   truncate-lines: t
   End:
