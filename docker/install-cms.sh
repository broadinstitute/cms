#!/bin/bash
#
# This script requires INSTALL_PATH (typically /opt/cms),
# CMS_PATH (typically /opt/cms/source), and
# CONDA_DEFAULT_ENV (typically /opt/miniconda) to be set.
#
# A miniconda install must exist at $CONDA_DEFAULT_ENV
# and $CONDA_DEFAULT_ENV/bin must be in the PATH
#
# Otherwise, this only requires the existence of the following files:
#	requirements-minimal.txt
#	requirements-conda.txt
#	requirements-conda-tests.txt
#	requirements-py3.txt

set -e -o pipefail -x

echo "PATH:              ${PATH}"
echo "INSTALL_PATH:      ${INSTALL_PATH}"
echo "CONDA_PREFIX:      ${CONDA_PREFIX}"
echo "CMS_PATH:    ${CMS_PATH}"
echo "MINICONDA_PATH:    ${MINICONDA_PATH}"
echo "CONDA_DEFAULT_ENV: ${CONDA_DEFAULT_ENV}"
CONDA_CHANNEL_STRING="--override-channels -c conda-forge -c bioconda -c defaults -c notestaff"

# setup/install cms directory tree and conda dependencies
sync

# manually install it ourselves instead of using easy-deploy
conda install -y \
      -q $CONDA_CHANNEL_STRING \
      --file "$CMS_PATH/requirements-cms.txt" \
      -p "${CONDA_PREFIX}"

# clean up
conda clean -y --all
