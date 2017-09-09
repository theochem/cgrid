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

    cdef cppclass Ln(ScalarFunction):
        Ln(double prefac, double alpha) except +
        double prefac()
        double alpha()

    cdef cppclass Linear(ScalarFunction):
        Linear(double slope, double offset) except +
        double slope()
        double offset()

    cdef cppclass Identity(ScalarFunction):
        Identity() except +

    cdef cppclass Constant(ScalarFunction):
        Constant(double offset) except +
        double offset()

    cdef cppclass Power(ScalarFunction):
        Power(double prefac, double power) except +
        double prefac()
        double power()

    cdef cppclass Rational(ScalarFunction):
        Rational(double prefac, double root) except +
        double prefac()
        double root()

    cdef cppclass Spline(ScalarFunction):
        Spline(size_t npoint) except +
        size_t npoint()

        double* values();
        double* derivs();
        double* derivs2();

        double left();
        double right();
        double x(const size_t index);

        void fit_derivs();

    cdef cppclass UniformCubicSpline(Spline):
        UniformCubicSpline(size_t npoint);
        UniformCubicSpline(size_t npoint, const double* const values);
        UniformCubicSpline(size_t npoint, const double* const values, double* const derivs);

    cdef cppclass Composed(ScalarFunction):
        Composed(Spline* spline, ScalarFunction* x_transform, ScalarFunction* y_transform,
                 ScalarFunction* left_extra, ScalarFunction* right_extra,
                 const double* values);

        Composed(Spline* spline, ScalarFunction* x_transform, ScalarFunction* y_transform,
                 ScalarFunction* left_extra, ScalarFunction* right_extra,
                 const double* values, const double* derivs);

        Spline* spline();
        ScalarFunction* x_transform();
        ScalarFunction* y_transform();
        ScalarFunction* left_extra();
        ScalarFunction* right_extra();
