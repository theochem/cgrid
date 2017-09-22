|Travis|
|Conda|
|PythonConda|
|Codecov|
|CondaVersion|

QCGrids is a numerical integration library for quantum chemistry.

Dependencies
============

-  celllists, python-celllists
-  C++11 compiler (tested: GNU and Intel)
-  Python, >=2.7.x, <3.x
-  Numpy, >1.9
-  Cython, >= 0.24.1

Runtime environment configuration
=================================

The instructions below explain how to install everything in
``${HOME}/.local``, such that you do not need root permissions to
install all the software. To make this work, some environment variables
must be set, e.g. in your ``~/.bashrc``.

::

    export PATH=${HOME}/.local/bin:${PATH}
    export LD_LIBRARY_PATH=${HOME}/.local/lib:${HOME}/.local/lib64:${LD_LIBRARY_PATH}

Installation of ``qcgrids``
===========================

Build (for installation in home directory):

::

    mkdir build
    cd build
    cmake .. -DCMAKE_BUILD_TYPE=release -DCMAKE_INSTALL_PREFIX=${HOME}/.local
    make

Install:

::

    make install

Installation of ``python-qcgrids``
==================================

In-place build and test

::

    cd python-qcgrids
    ./setup.py build_ext -i -L${LD_LIBRARY_PATH}
    nosetests -v

Build and install (into home directory):

::

    ./setup.py build_ext -L${LD_LIBRARY_PATH}
    ./setup.py install --user
    
.. |Travis| image:: https://travis-ci.org/theochem/qcgrids.svg?branch=master
    :target: https://travis-ci.org/theochem/qcgrids
.. |Codecov| image:: https://img.shields.io/codecov/c/github/theochem/qcgrids/master.svg
    :target: https://codecov.io/gh/theochem/qcgrids
.. |Conda| image:: https://img.shields.io/conda/v/theochem/qcgrids.svg
    :target: https://anaconda.org/theochem/qcgrids
.. |PythonConda| image:: https://img.shields.io/conda/vn/theochem/python-qcgrids.svg
    :target: https://anaconda.org/theochem/python-qcgrids
.. |CondaVersion| image:: https://img.shields.io/conda/pn/theochem/qcgrids.svg
    :target: https://anaconda.org/theochem/qcgrids

