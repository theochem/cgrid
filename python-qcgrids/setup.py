#!/usr/bin/env python
# -*- coding: utf-8 -*-
# QCGrids is a numerical integration library for quantum chemistry.
# Copyright (C) 2011-2017 The QCGrids developers
#
# This file is part of QCGrids.
#
# QCGrids is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# QCGrids is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>
# --

import os
import numpy as np
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext


def get_version():
    """Load the version from version.py, without importing it.

    This function assumes that the last line in the file contains a variable defining the
    version string with single quotes.

    """
    with open('qcgrids/version.py', 'r') as f:
        return f.read().split('=')[-1].replace('\'', '').strip()


def get_cxxflags():
    """If the CXXFLAGS variable is defined (clang/osx) then get it."""
    return os.environ.get("CXXFLAGS", "").split()


def get_include_path():
    """If the PREFIX variable is defined (conda) then get the conda include prefix."""
    prefix = os.environ.get("PREFIX", "")
    if prefix:
        return [os.path.join(prefix, "include")]
    else:
        return []


setup(
    name='python-qcgrids',
    version=get_version(),
    description='QCGrids is a numerical integration library for quantum chemistry.',
    author='Toon Verstraelen',
    author_email='Toon.Verstraelen@UGent.be',
    cmdclass = {'build_ext': build_ext},
    packages = ['qcgrids'],
    package_data = {
        'qcgrids': ['ext.pxd', 'cellgrid.pxd', 'scalarfns.pxd'],
    },
    ext_modules=[
        Extension(
            "qcgrids.ext",
            sources=['qcgrids/ext.pyx'],
            depends=['qcgrids/ext.pxd', 'qcgrids/cellgrid.pxd', 'qcgrids/scalarfns.pxd'],
            libraries=['qcgrids'],
            include_dirs=[np.get_include()] + get_include_path(),
            extra_compile_args=get_cxxflags() or ['-std=c++11', '-Wall'],
            language="c++")],
)
