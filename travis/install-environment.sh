#!/bin/bash

if [[ "$TRAVIS_PYTHON_VERSION" == 2* ]]; then
    conda env create -f conda-environment_py2.yml -p $HOME/cms-env
else
    conda env create -f conda-environment_py3.yml -p $HOME/cms-env
fi

