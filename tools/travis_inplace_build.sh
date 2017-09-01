#!/usr/bin/env bash
set -ev

# In-place building and testing
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This is allows for coverage analysis and dynamic linting. The compiler settings used
# here are not suitable for releases, so we need to recompile and rerun the tests.

# Install GCC compilers for in-place builds, even on OSX because clang does not manage to
# compiler our C++ code.
conda install gcc libgcc;

# PY project
# Uninstall conda package, to be sure. The conda cpp package is still used.
conda uninstall ${CONDA_PKG_NAME_PY}
(cd ${PYDIR}; python setup.py build_ext -i --define CYTHON_TRACE_NOGIL)
# Run nosetests without coverage.xml output. That file is broken by nosetests (pyx files
# not include) and gets priority over .coverage, which contains everything.
(cd ${PYDIR}; nosetests ${PROJECT_NAME} -v --detailed-errors --with-coverage --cover-package=${PROJECT_NAME} --cover-tests --cover-erase --cover-inclusive --cover-branches; coverage xml -i)

if [ "$TRAVIS_PULL_REQUEST" != "false" ]; then
  (cd ${PYDIR}; cardboardlinter --refspec $TRAVIS_BRANCH -f 'dynamic');
fi

# CPP project
# Uninstall conda package, to be sure.
conda uninstall ${CONDA_PKG_NAME_CPP}

# Install dependencies
conda install ${DEPENDENCIES}
(mkdir debug; cd debug; cmake cmake -DCMAKE_BUILD_TYPE=debug ..; make)
# DYLD_LIBRARY_PATH is needed on OSX to let it find the shared libraries in
# ${CONDA_PREFIX}/lib. It should not have any effect on Linux.
(cd debug; DYLD_LIBRARY_PATH=${DYLD_LIBRARY_PATH}:${CONDA_PREFIX}/lib/ ./${PROJECT_NAME}/tests/test_${PROJECT_NAME})

if [ "$TRAVIS_PULL_REQUEST" != "false" ]; then
  cardboardlinter --refspec $TRAVIS_BRANCH -f 'dynamic';
fi