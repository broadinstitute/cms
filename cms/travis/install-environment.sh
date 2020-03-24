#!/bin/bash

if [[ "$(python -c 'import sys; print(sys.version[0])')" == 2* ]]; then
    echo "Creating conda environment using Py2 environment file"
    conda env create -f conda-environment_py2.yaml -p $HOME/cms-env
else
    echo "Creating conda environment using Py3 environment file"
    conda env create -f conda-environment_py3.yaml -p $HOME/cms-env
fi

