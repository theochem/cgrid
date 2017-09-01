#!/usr/bin/env bash
set -ev

# Conda-based tests
# ~~~~~~~~~~~~~~~~~
# This is needed to simulate the end-user situation, and it also creates packages for a
# conda release.

# Build the conda packages
conda build -q tools/conda.recipe
(cd ${PYDIR}; conda build -q tools/conda.recipe)
# Install Conda packages
conda install --use-local ${CONDA_PKG_NAME_CPP} ${CONDA_PKG_NAME_PY}
# Run the unittests, using the installed package
${CONDA_PREFIX}/libexec/${PROJECT_NAME}/test_${PROJECT_NAME}
(cd; nosetests ${PROJECT_NAME} -v --detailed-errors)
