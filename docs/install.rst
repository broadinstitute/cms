Installation
============


System dependencies
-------------------

To be described in greater detail...

Manual Installation
-------------------


Step 1: Install Conda

To use conda, you need to install the `Conda package manager <http://conda.pydata.org/miniconda.html>`_ which is most easily obtained via the Miniconda Python distribution. Miniconda can be installed to your home directory without admin priviledges. On Broad Institute systems, you can make use of the ".anaconda3-4.0.0" dotkit.

Step 2: Configure Conda

Software used by the cms project is distributed through the bioconda channel for the conda package manager. It is necessary to add this channel to the conda config::

  conda config --add channels bioconda

Step 3: Make a conda environment and install cms

It is recommended to install cms into its own conda directory. This ensures its dependencies do not interfere with other conda packages installed on your system. A new conda environment can be created with the following command, which will also install relevant cms dependencies::

  conda env create -f=conda-environment.yml -n cms-env

Step 4: Activate the cms environment

In order to use cms, you will need to activate its conda environment::

  source activate cms-env

