.. image:: https://travis-ci.org/theochem/cgrid.svg?branch=master
    :target: https://travis-ci.org/theochem/cgrid
.. image:: https://img.shields.io/codecov/c/github/theochem/cgrid/master.svg
    :target: https://codecov.io/gh/theochem/cgrid
.. image:: https://img.shields.io/conda/v/theochem/cgrid.svg
    :target: https://anaconda.org/theochem/cgrid
.. image:: https://img.shields.io/conda/vn/theochem/python-cgrid.svg
    :target: https://anaconda.org/theochem/python-cgrid
.. image:: https://img.shields.io/conda/pn/theochem/cgrid.svg
    :target: https://anaconda.org/theochem/cgrid
.. image:: https://img.shields.io/github/release/theochem/cgrid.svg
    :target: https://github.com/theochem/cgrid/releases

CGrid is a library for for molecular numerical integration.

Installation
============

When you are interested in using cgrid (without needing to modify it), you
can install cgrid with conda. After installing and activating a miniconda
environment, run:

.. code-block:: bash

  conda install -c theochem cgrid python-cgrid

If you are interesed in working on the development of cgrid, you first need
to check out the latest version from the git repository

.. code-block:: bash

  git clone git@github.com:theochem/cgrid.git
  cd cgrid

Then install Roberto and run it in the root of the repository:

.. code-block:: bash

  pip install --user --upgrade 'roberto<2.0.0'
  rob quality

This will build cgrid in-place and run all tests. More details for
potential contributors are given in `CONTRIBUTING.rst <CONTRIBUTING.rst>`_.
