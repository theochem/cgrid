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
# cython: linetrace=True, embedsignature=True, language_level=3


from libcpp.memory cimport unique_ptr

cimport cgrid.cellgrid as cellgrid
cimport cgrid.scalarfns as scalarfns


cdef class GridFunc:
    cdef cellgrid.GridFunc _funcptr
    cdef void* _extra_arg


cdef class Cellgrid:
    cdef cellgrid.Cellgrid* _this


cdef class ScalarFunction:
    cdef unique_ptr[scalarfns.ScalarFunction] _this


cdef class Exp(ScalarFunction):
    pass


cdef class Ln(ScalarFunction):
    pass


cdef class Linear(ScalarFunction):
    pass


cdef class Identity(ScalarFunction):
    pass


cdef class Constant(ScalarFunction):
    pass


cdef class Power(ScalarFunction):
    pass


cdef class Rational(ScalarFunction):
    pass


cdef class Spline(ScalarFunction):
    pass


cdef class UniformCubicSpline(Spline):
    pass


cdef class Composed(ScalarFunction):
    cdef Spline _spline
    cdef ScalarFunction _x_transform
    cdef ScalarFunction _y_transform
    cdef ScalarFunction _left_extra
    cdef ScalarFunction _right_extra
