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


from libcpp.memory cimport unique_ptr

cimport cellgrid
cimport scalarfns


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
