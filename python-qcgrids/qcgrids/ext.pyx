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
# cython: linetrace=True
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
    'ScalarFunction', 'Exp',
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
    def __init__(self):
        raise RuntimeError("Cannot instantiate abstract base class ScalarFunction.")

    property invertible:
        def __get__(self):
            return deref(self._this).invertible()

    def calc(self, xpy, int nderiv):
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
        result = np.zeros_like(x)
        cdef np.broadcast it = np.broadcast(x, result)
        while np.PyArray_MultiIter_NOTDONE(it):
            (<double*>np.PyArray_MultiIter_DATA(it, 1))[0] = \
                deref(self._this).value(deref(<double*>np.PyArray_MultiIter_DATA(it, 0)))
            np.PyArray_MultiIter_NEXT(it)
        return result if isinstance(x, np.ndarray) else result[()]

    def value_inv(self, x):
        result = np.zeros_like(x)
        cdef np.broadcast it = np.broadcast(x, result)
        while np.PyArray_MultiIter_NOTDONE(it):
            (<double*>np.PyArray_MultiIter_DATA(it, 1))[0] = \
                deref(self._this).value_inv(deref(<double*>np.PyArray_MultiIter_DATA(it, 0)))
            np.PyArray_MultiIter_NEXT(it)
        return result if isinstance(x, np.ndarray) else result[()]

    def deriv(self, x):
        result = np.zeros_like(x)
        cdef np.broadcast it = np.broadcast(x, result)
        while np.PyArray_MultiIter_NOTDONE(it):
            (<double*>np.PyArray_MultiIter_DATA(it, 1))[0] = \
                deref(self._this).deriv(deref(<double*>np.PyArray_MultiIter_DATA(it, 0)))
            np.PyArray_MultiIter_NEXT(it)
        return result if isinstance(x, np.ndarray) else result[()]

    def deriv_inv(self, x):
        result = np.zeros_like(x)
        cdef np.broadcast it = np.broadcast(x, result)
        while np.PyArray_MultiIter_NOTDONE(it):
            (<double*>np.PyArray_MultiIter_DATA(it, 1))[0] = \
                deref(self._this).deriv_inv(deref(<double*>np.PyArray_MultiIter_DATA(it, 0)))
            np.PyArray_MultiIter_NEXT(it)
        return result if isinstance(x, np.ndarray) else result[()]

    def deriv2(self, x):
        result = np.zeros_like(x)
        cdef np.broadcast it = np.broadcast(x, result)
        while np.PyArray_MultiIter_NOTDONE(it):
            (<double*>np.PyArray_MultiIter_DATA(it, 1))[0] = \
                deref(self._this).deriv2(deref(<double*>np.PyArray_MultiIter_DATA(it, 0)))
            np.PyArray_MultiIter_NEXT(it)
        return result if isinstance(x, np.ndarray) else result[()]

    def deriv2_inv(self, x):
        result = np.zeros_like(x)
        cdef np.broadcast it = np.broadcast(x, result)
        while np.PyArray_MultiIter_NOTDONE(it):
            (<double*>np.PyArray_MultiIter_DATA(it, 1))[0] = \
                deref(self._this).deriv2_inv(deref(<double*>np.PyArray_MultiIter_DATA(it, 0)))
            np.PyArray_MultiIter_NEXT(it)
        return result if isinstance(x, np.ndarray) else result[()]


cdef class Exp(ScalarFunction):
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
