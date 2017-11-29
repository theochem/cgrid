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
#
# Needed for coverage analysis:
# cython: linetrace=True, embedsignature=True
'''Python wrapper for the QCGrids library'''


from cpython.ref cimport PyTypeObject
from libc.string cimport memcpy
from cpython cimport Py_INCREF
from cython.operator cimport dereference as deref


import numpy as np
cimport numpy as np
np.import_array()

cimport cellcutoff.ext as cellcutoff

cimport cellgrid
cimport scalarfns


cdef extern from "numpy/arrayobject.h":
    object PyArray_NewFromDescr(PyTypeObject* subtype, np.dtype descr,
                                int nd, np.npy_intp* dims, np.npy_intp* strides,
                                void* data, int flags, object obj)


__all__ = [
    # cellgrid.h
    'Cellgrid'
    # scalarfns.h
    'ScalarFunction', 'Exp', 'Ln', 'Linear', 'Identity', 'Constant', 'Power', 'Rational',
    'Spline', 'UniformCubicSpline', 'Composed',
]


def check_array_arg(name, arg, expected_shape):
    if not arg.flags['C_CONTIGUOUS']:
        raise TypeError('Argument %s must be C_CONTIHUOUS.' % arg)
    for i, n in enumerate(expected_shape):
        if n >= 0 and arg.shape[i] != n:
            raise TypeError(('Axis %i of argument %s has length %i but while '
                             'expecting %i') % (i, name, arg.shape[i], n))


cdef class Cellgrid(object):
    def __cinit__(self, cellcutoff.Cell cell, double spacing=1):
        self._this = new cellgrid.Cellgrid(cell._this[0], spacing)

    def __init__(self, cellcutoff.Cell cell, double spacing=1):
        pass

    def __dealloc__(self):
        if self._this != NULL:
            del self._this

    property npoint:
        def __get__(self):
            return self._this.grid_array().size()

    property cell:
        def __get__(self):
            cdef np.ndarray[double, ndim=2] vecs = np.zeros((3, 3), float)
            memcpy(&vecs[0, 0], self._this.cell().vecs(), sizeof(double)*9);
            return cellcutoff.Cell(vecs)

    property subcell:
        def __get__(self):
            cdef np.ndarray[double, ndim=2] vecs = np.zeros((3, 3), float)
            memcpy(&vecs[0, 0], self._this.subcell().vecs(), sizeof(double)*9);
            return cellcutoff.Cell(vecs)

    property shape:
        def __get__(self):
            cdef np.ndarray[long, ndim=1] shape = np.zeros(3, int)
            memcpy(&shape[0], self._this.shape(), sizeof(int)*3);
            return shape

    property grid_array:
        def __get__(self):
            cdef np.npy_intp dims[1]
            dims[0] = self._this.grid_array().size()
            cdef np.dtype dtype = np.dtype([
                ('cart', np.double, 3),
                ('icell', np.intc, 3),
                ('weight', np.double),
                ('index', np.intc),
            ], align=True)
            assert dtype.itemsize == sizeof(cellgrid.CellgridPoint)
            Py_INCREF(dtype)
            result = PyArray_NewFromDescr(
                <PyTypeObject*> np.ndarray,
                dtype,
                1,
                dims,
                NULL,
                <void*> self._this.grid_array().data(),
                0,
                None)
            np.set_array_base(result, self)
            return result

    property points:
        def __get__(self):
            return self.grid_array['cart']

    property weights:
        def __get__(self):
            return self.grid_array['weight']

    property indices:
        def __get__(self):
            return self.grid_array['index']

    def __len__(self):
        return self._this.grid_array().size()

    def append_many(self, np.ndarray[double, ndim=2] points not None,
                    np.ndarray[double, ndim=1] weights not None):
        cdef size_t npoint = 0;
        check_array_arg('points', points, (-1, 3))
        check_array_arg('weights', weights, (-1,))
        if points.shape[0] != weights.shape[0]:
            raise TypeError('Points and weights must have the same length')
        npoint = points.shape[0]
        self._this.emplace_back_many(&points[0, 0], &weights[0], npoint)

    def sort(self):
        self._this.sort()

    def iadd_cutoff(self, np.ndarray[double, ndim=1] center not None,
                    double cutoff, GridFunc grid_func not None,
                    np.ndarray[double, ndim=1] output not None):
        check_array_arg('center', center, (3,))
        check_array_arg('output', output, (-1,))
        self._this.iadd_cutoff(&center[0], cutoff, grid_func._funcptr,
                               grid_func._extra_arg, &output[0])

    def integrate_cutoff(self, np.ndarray[double, ndim=1] center not None,
                         double cutoff, GridFunc grid_func not None,
                         np.ndarray[double, ndim=1] factor not None):
        check_array_arg('center', center, (3,))
        check_array_arg('factor', factor, (-1,))
        return self._this.integrate_cutoff(&center[0], cutoff, grid_func._funcptr,
                                           grid_func._extra_arg, &factor[0])

    def integrate(self, *factors):
        subscripts = '%s->' % (','.join(['i']*(len(factors)+1)))
        return np.einsum(subscripts, self.weights, *factors)


cdef class ScalarFunction(object):
    """Abstract base class for all scalar functions."""

    def __init__(self):
        raise RuntimeError("Cannot instantiate abstract base class ScalarFunction.")

    property invertible:
        def __get__(self):
            """Is the function invertible or not?"""
            return deref(self._this).invertible()

    def calc(self, xpy, int nderiv):
        """Compute the function value and optional first and second derivatives.

        Parameters
        ----------
        xpy : float or np.ndarray
            The input to the function.
        nderiv : int
            The number of derivatives (0, 1, 2).

        Returns
        -------
        output : np.ndarray
            Array with ``shape = xpy.shape + (nderiv + 1,)``. The last index is used for
            function value (0), first derivative (1) and second derivative (3).
        """
        x = np.asarray(xpy, dtype=float)
        output = np.zeros(x.shape + (nderiv + 1,), dtype=float)
        cdef np.flatiter itx = np.PyArray_IterNew(x)
        cdef int xndim = x.ndim
        cdef np.flatiter ito = np.PyArray_IterAllButAxis(output, &(xndim))
        while np.PyArray_ITER_NOTDONE(itx):
            deref(self._this).calc(
                deref(<double*>np.PyArray_ITER_DATA(itx)),
                nderiv, <double*>np.PyArray_ITER_DATA(ito),
            )
            np.PyArray_ITER_NEXT(itx)
            np.PyArray_ITER_NEXT(ito)
        return output

    def calc_inv(self, xpy, int nderiv):
        """Compute the inverse function value and optional first and second derivatives.

        Parameters
        ----------
        xpy : float or np.ndarray
            The input to the function.
        nderiv : int
            The number of derivatives (0, 1, 2).

        Returns
        -------
        output : np.ndarray
            Array with ``shape = xpy.shape + (nderiv + 1,)``. The last index is used for
            function value (0), first derivative (1) and second derivative (3).
        """
        x = np.asarray(xpy, dtype=float)
        output = np.zeros(x.shape + (nderiv + 1,), dtype=float)
        cdef np.flatiter itx = np.PyArray_IterNew(x)
        cdef int xndim = x.ndim
        cdef np.flatiter ito = np.PyArray_IterAllButAxis(output, &(xndim))
        while np.PyArray_ITER_NOTDONE(itx):
            deref(self._this).calc_inv(
                deref(<double*>np.PyArray_ITER_DATA(itx)),
                nderiv, <double*>np.PyArray_ITER_DATA(ito),
            )
            np.PyArray_ITER_NEXT(itx)
            np.PyArray_ITER_NEXT(ito)
        return output

    def value(self, x):
        """Compute the function value. Works as u-func."""
        result = np.zeros_like(x, dtype=float)
        cdef np.broadcast it = np.broadcast(np.asarray(x, dtype=float), result)
        while np.PyArray_MultiIter_NOTDONE(it):
            (<double*>np.PyArray_MultiIter_DATA(it, 1))[0] = \
                deref(self._this).value(deref(<double*>np.PyArray_MultiIter_DATA(it, 0)))
            np.PyArray_MultiIter_NEXT(it)
        return result if isinstance(x, np.ndarray) else result[()]

    def value_inv(self, x):
        """Compute the inverse function value. Works as u-func."""
        result = np.zeros_like(x, dtype=float)
        cdef np.broadcast it = np.broadcast(np.asarray(x, dtype=float), result)
        while np.PyArray_MultiIter_NOTDONE(it):
            (<double*>np.PyArray_MultiIter_DATA(it, 1))[0] = \
                deref(self._this).value_inv(deref(<double*>np.PyArray_MultiIter_DATA(it, 0)))
            np.PyArray_MultiIter_NEXT(it)
        return result if isinstance(x, np.ndarray) else result[()]

    def deriv(self, x):
        """Compute the derivative. Works as u-func."""
        result = np.zeros_like(x, dtype=float)
        cdef np.broadcast it = np.broadcast(np.asarray(x, dtype=float), result)
        while np.PyArray_MultiIter_NOTDONE(it):
            (<double*>np.PyArray_MultiIter_DATA(it, 1))[0] = \
                deref(self._this).deriv(deref(<double*>np.PyArray_MultiIter_DATA(it, 0)))
            np.PyArray_MultiIter_NEXT(it)
        return result if isinstance(x, np.ndarray) else result[()]

    def deriv_inv(self, x):
        """Compute the derivative of the inverse. Works as u-func."""
        result = np.zeros_like(x, dtype=float)
        cdef np.broadcast it = np.broadcast(np.asarray(x, dtype=float), result)
        while np.PyArray_MultiIter_NOTDONE(it):
            (<double*>np.PyArray_MultiIter_DATA(it, 1))[0] = \
                deref(self._this).deriv_inv(deref(<double*>np.PyArray_MultiIter_DATA(it, 0)))
            np.PyArray_MultiIter_NEXT(it)
        return result if isinstance(x, np.ndarray) else result[()]

    def deriv2(self, x):
        """Compute the second derivative. Works as u-func."""
        result = np.zeros_like(x, dtype=float)
        cdef np.broadcast it = np.broadcast(np.asarray(x, dtype=float), result)
        while np.PyArray_MultiIter_NOTDONE(it):
            (<double*>np.PyArray_MultiIter_DATA(it, 1))[0] = \
                deref(self._this).deriv2(deref(<double*>np.PyArray_MultiIter_DATA(it, 0)))
            np.PyArray_MultiIter_NEXT(it)
        return result if isinstance(x, np.ndarray) else result[()]

    def deriv2_inv(self, x):
        """Compute the second derivative of the inverse. Works as u-func."""
        result = np.zeros_like(x, dtype=float)
        cdef np.broadcast it = np.broadcast(np.asarray(x, dtype=float), result)
        while np.PyArray_MultiIter_NOTDONE(it):
            (<double*>np.PyArray_MultiIter_DATA(it, 1))[0] = \
                deref(self._this).deriv2_inv(deref(<double*>np.PyArray_MultiIter_DATA(it, 0)))
            np.PyArray_MultiIter_NEXT(it)
        return result if isinstance(x, np.ndarray) else result[()]


cdef class Exp(ScalarFunction):
    """Exponential function: y = exp(slope*x + offset)."""

    def __cinit__(self, double slope, double offset):
        self._this.reset(new scalarfns.Exp(slope, offset))

    def __init__(self, double slope, double offset):
        pass

    property slope:
        def __get__(self):
            return deref(<scalarfns.Exp*?>self._this.get()).slope()

    property offset:
        def __get__(self):
            return deref(<scalarfns.Exp*?>self._this.get()).offset()


cdef class Ln(ScalarFunction):
    """Logarithmic function: y = prefac*ln(alpha*x)."""

    def __cinit__(self, double prefac, double alpha):
        self._this.reset(new scalarfns.Ln(prefac, alpha))

    def __init__(self, double prefac, double alpha):
        pass

    property prefac:
        def __get__(self):
            return deref(<scalarfns.Ln*?>self._this.get()).prefac()

    property alpha:
        def __get__(self):
            return deref(<scalarfns.Ln*?>self._this.get()).alpha()


cdef class Linear(ScalarFunction):
    """Linear function: y = slone*x + offset."""

    def __cinit__(self, double slope, double offset):
        self._this.reset(new scalarfns.Linear(slope, offset))

    def __init__(self, double slope, double offset):
        pass

    property slope:
        def __get__(self):
            return deref(<scalarfns.Linear*?>self._this.get()).slope()

    property offset:
        def __get__(self):
            return deref(<scalarfns.Linear*?>self._this.get()).offset()


cdef class Identity(ScalarFunction):
    """Identity function: y = x."""

    def __cinit__(self):
        self._this.reset(new scalarfns.Identity())

    def __init__(self):
        pass


cdef class Constant(ScalarFunction):
    """Constant function: y = offset."""

    def __cinit__(self, double offset):
        self._this.reset(new scalarfns.Constant(offset))

    def __init__(self, double offset):
        pass

    property offset:
        def __get__(self):
            return deref(<scalarfns.Constant*?>self._this.get()).offset()


cdef class Power(ScalarFunction):
    """Power function: y = prefac*x**power."""

    def __cinit__(self, double prefac, double power):
        self._this.reset(new scalarfns.Power(prefac, power))

    def __init__(self, double prefac, double power):
        pass

    property prefac:
        def __get__(self):
            return deref(<scalarfns.Power*?>self._this.get()).prefac()

    property power:
        def __get__(self):
            return deref(<scalarfns.Power*?>self._this.get()).power()


cdef class Rational(ScalarFunction):
    """Rational function: y = prefac*x/(1 + x/root)."""

    def __cinit__(self, double prefac, double root):
        self._this.reset(new scalarfns.Rational(prefac, root))

    def __init__(self, double prefac, double root):
        pass

    property prefac:
        def __get__(self):
            return deref(<scalarfns.Rational*?>self._this.get()).prefac()

    property root:
        def __get__(self):
            return deref(<scalarfns.Rational*?>self._this.get()).root()


cdef class Spline(ScalarFunction):
    """Abstract base class for all splines."""

    def __init__(self, npoint):
        raise RuntimeError("Cannot instantiate abstract base class Spline.")

    property npoint:
        def __get__(self):
            """Number of grid points."""
            return deref(<scalarfns.Spline*?>self._this.get()).npoint()

    property values:
        def __get__(self):
            """Function values at the grid points that parameterize the spline."""
            cdef double* ptr = deref(<scalarfns.Spline*?>self._this.get()).values()
            if ptr != NULL:
                return np.asarray(<np.float64_t[:self.npoint:1]>ptr)

    property derivs:
        def __get__(self):
            """Derivatives values at the grid points that parameterize the spline."""
            cdef double* ptr = deref(<scalarfns.Spline*?>self._this.get()).derivs()
            if ptr != NULL:
                return np.asarray(<np.float64_t[:self.npoint:1]>ptr)

    property derivs2:
        def __get__(self):
            """Second derivatives values at the grid points that parameterize the spline."""
            cdef double* ptr = deref(<scalarfns.Spline*?>self._this.get()).derivs2()
            if ptr != NULL:
                return np.asarray(<np.float64_t[:self.npoint:1]>ptr)

    property left:
        def __get__(self):
            """Left-most grid point."""
            return deref(<scalarfns.Spline*?>self._this.get()).left()

    property right:
        def __get__(self):
            """Right-most grid point."""
            return deref(<scalarfns.Spline*?>self._this.get()).right()

    def x(self, size_t index):
        """Position of grid point with given index, in [0, npoint - 1]."""
        return deref(<scalarfns.Spline*?>self._this.get()).x(index)

    def fit_derivs(self):
        """Fit the derivatives of a natural cubic spline, using the function values."""
        deref(<scalarfns.Spline*?>self._this.get()).fit_derivs()


cdef class UniformCubicSpline(Spline):
    """Cubic spline with grid points [0.0, 1.0, ... npoint - 1]."""

    def __cinit__(self, np.float64_t[::1] values not None, np.float64_t[::1] derivs=None):
        cdef size_t npoint = values.shape[0]
        if derivs is None:
            self._this.reset(new scalarfns.UniformCubicSpline(npoint, &values[0]))
        else:
            if derivs.shape[0] != npoint:
                raise TypeError("The arrays values and derivs must have the same size.")
            self._this.reset(new scalarfns.UniformCubicSpline(npoint, &values[0], &derivs[0]))

    def __init__(self, np.float64_t[::1] values=None, np.float64_t[::1] derivs=None):
        pass


cdef class Composed(ScalarFunction):
    """Composition of spline and additional functions.

    Central to the composed function is a spline, for now only an instance
    UniformCubicSpline. For a given value x, the composed function is computed as:

    - if x < x_transform.value_inv(0) : right_extra(x)
    - if x > x_transform.value_inv(npoint - 1) : left_extra(x)
    - else : y_transform.value(spline.value(x_transform.value_inv(x)))
    """

    def __cinit__(self, ScalarFunction x_transform, ScalarFunction y_transform,
                  ScalarFunction left_extra, ScalarFunction right_extra,
                  np.float64_t[::1] values not None, np.float64_t[::1] derivs=None):
        cdef size_t npoint = values.shape[0]
        self._spline = UniformCubicSpline(np.zeros(npoint, dtype=float))
        self._x_transform = x_transform
        self._y_transform = y_transform
        self._left_extra = left_extra
        self._right_extra = right_extra
        if derivs is None:
            self._this.reset(new scalarfns.Composed(
                <scalarfns.Spline*?>self._spline._this.get(),
                x_transform._this.get() if x_transform is not None else NULL,
                y_transform._this.get() if y_transform is not None else NULL,
                left_extra._this.get() if left_extra is not None else NULL,
                right_extra._this.get() if right_extra is not None else NULL,
                &values[0]))
        else:
            if derivs.shape[0] != npoint:
                raise TypeError("The arrays values and derivs must have the same size.")
            self._this.reset(new scalarfns.Composed(
                <scalarfns.Spline*?>self._spline._this.get(),
                x_transform._this.get() if x_transform is not None else NULL,
                y_transform._this.get() if y_transform is not None else NULL,
                left_extra._this.get() if left_extra is not None else NULL,
                right_extra._this.get() if right_extra is not None else NULL,
                &values[0], &derivs[0]))

    def __init__(self, ScalarFunction x_transform, ScalarFunction y_transform,
                  ScalarFunction left_extra, ScalarFunction right_extra,
                  np.float64_t[::1] values not None, np.float64_t[::1] derivs=None):
        pass

    property spline:
        def __get__(self):
            """Spline used in the composed function."""
            return self._spline

    property x_transform:
        def __get__(self):
            """Transformation from spline x coordinate to composed x coordinate."""
            return self._x_transform

    property y_transform:
        def __get__(self):
            """Transformation from spline y coordinate to composed y coordinate."""
            return self._y_transform

    property left_extra:
        def __get__(self):
            """Extrapolation to the left of the spline."""
            return self._left_extra

    property right_extra:
        def __get__(self):
            """Extrapolation to the right of the spline."""
            return self._right_extra
