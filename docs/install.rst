Installation
============


System dependencies
-------------------

To be described in greater detail...

Python dependencies
-------------------

The **command line tools** require Python >= 2.7. Required packages are listed in requirements.txt and can be
installed the usual pip way::

  pip install -r requirements.txt

Additionally, in order to use the **pipeline infrastructure**, Python 3.4
is required (Python 2 is not supported) and you must install snakemake
as well::

  pip install -r requirements-pipes.txt

However, most of the real functionality is encapsulated in the command line
tools, which can be used without any of the pipeline infrastructure.

You should either sudo pip install or use a virtualenv (recommended).


Tool dependencies
-----------------

A lot of effort has gone into writing auto download/compile wrappers for
most of the bioinformatic tools we rely on here. They will auto-download
and install the first time they are needed by any command. If you want
to pre-install all of the external tools, simply type this::

  python -m unittest test.test_tools.TestToolsInstallation -v
