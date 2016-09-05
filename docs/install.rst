Installation
============

Step 1: Install Conda

To use conda, you need to install the `Conda package manager <http://conda.pydata.org/miniconda.html>`_ which is most easily obtained via the Miniconda Python distribution. Miniconda can be installed to your home directory without admin priviledges. 

Step 2: Configure Conda

Software used by the cms project is distributed through the bioconda channel for the conda package manager. It is necessary to add this channel to the conda config::

  conda config --add channels bioconda

Step 3: Make a conda environment and install cms

It is recommended to install cms into its own conda directory. This ensures its dependencies do not interfere with other conda packages installed on your system. First clone the source repository from Github::
	git clone git@github.com:broadinstitute/cms.git

A new conda environment can be created with the following command, which will also install relevant cms dependencies. It is recommended to use the Python3 version of the environment file::

  conda env create -f=conda-environment_py3.yaml -n cms-env3

Step 4: Activate the cms environment

In order to use cms, you will need to activate its conda environment::

  source activate cms-env3

