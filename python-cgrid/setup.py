#!/usr/bin/env python3
# CGrid is a library for for molecular numerical integration.
# Copyright (C) 2011-2019 The CGrid Development Team
#
# This file is part of CGrid.
#
# CGrid is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# CGrid is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>
# --
"""Package build and install script."""

import os
import sys

from Cython.Distutils.build_ext import new_build_ext as build_ext
import numpy as np
from setuptools import setup, Extension


def get_version():
    """Read __version__ from version.py, with exec to avoid importing it."""
    try:
        with open(os.path.join('cgrid', 'version.py'), 'r') as f:
            myglobals = {}
            exec(f.read(), myglobals)  # pylint: disable=exec-used
        return myglobals['__version__']
    except IOError:
        return "0.0.0.post0"


# Because we follow distutils' naming conventions...
# pylint: disable=invalid-name
class my_build_ext(build_ext):
    """Workaround for rpath bug in distutils for OSX."""

    def finalize_options(self):
        super().finalize_options()
        # Special treatment of rpath in case of OSX, to work around python
        # distutils bug 36353. This constructs proper rpath arguments for clang.
        # See https://bugs.python.org/issue36353
        if sys.platform[:6] == "darwin":
            for path in self.rpath:
                for ext in self.extensions:
                    ext.extra_link_args.append("-Wl,-rpath," + path)
            self.rpath[:] = []


setup(
    name='python-cgrid',
    version=get_version(),
    description='CGrid is a library for for molecular numerical integration',
    author='The CGrid development team',
    cmdclass={'build_ext': build_ext},
    packages=['cgrid'],
    package_data={
        'cgrid': ['ext.pxd', 'cellgrid.pxd', 'scalarfns.pxd'],
    },
    ext_modules=[Extension(
        "cgrid.ext",
        sources=['cgrid/ext.pyx'],
        depends=['cgrid/ext.pxd', 'cgrid/cellgrid.pxd', 'cgrid/scalarfns.pxd'],
        libraries=['cgrid'],
        include_dirs=[np.get_include()],
        extra_compile_args=['-std=c++11', '-Wall'],
        language="c++"
    )],
    classifiers=[
        'Environment :: Console',
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Topic :: Scientific/Engineering :: Physics',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Intended Audience :: Science/Research',
    ],
    setup_requires=['numpy>=1.0', 'cython>=0.28.0'],
    install_requires=['numpy>=1.0', 'cython>=0.28.0', 'python-cellcutoff'],
)
