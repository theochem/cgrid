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


from libcpp cimport bool

cdef extern from "qcgrids/scalarfns.h" namespace "qcgrids":
    cdef cppclass ScalarFunction:
        bool invertible()
        void calc(const double x, const int nderiv, double* const output) except +
        void calc_inv(const double y, const int nderiv, double* const output) except +
        double value(const double x)
        double value_inv(const double y)
        double deriv(const double x)
        double deriv_inv(const double y)
        double deriv2(const double x)
        double deriv2_inv(const double y)

    cdef cppclass Exp(ScalarFunction):
        Exp(double slope, double offset) except +
        double slope()
        double offset()
