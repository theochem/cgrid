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


from nose.tools import assert_raises
import numpy as np
from numpy.testing import assert_allclose

from cellcutoff import Cell
from qcgrids.ext import Cellgrid, Exp


def test_something():
    cell = Cell(np.identity(3)*5.0)
    grid = Cellgrid(cell)
    grid.append_many(np.random.uniform(0, 10, (2, 3)), np.random.uniform(0, 1, 2))
    assert grid.npoint == 2


def test_exp():
    expfn = Exp(0.2, -0.3)
    assert expfn.slope == 0.2
    assert expfn.offset == -0.3
    check_fn(expfn, 1.2, np.exp(0.2*1.2 - 0.3))


def check_fn(func, x, y):
    """Test if a subclass of ScalarFunction has the right API from the base class."""
    # Test calc
    output = func.calc(x, 2)
    assert output.shape == (3,)
    assert_allclose(output[0], func.value(x))
    assert_allclose(output[1], func.deriv(x))
    assert_allclose(output[2], func.deriv2(x))
    output = func.calc(x, 1)
    assert output.shape == (2,)
    assert_allclose(output[0], func.value(x))
    assert_allclose(output[1], func.deriv(x))
    output = func.calc(x, 0)
    assert output.shape == (1,)
    assert_allclose(output[0], func.value(x))
    with assert_raises(ValueError):
        func.calc(x, 3)
    with assert_raises(ValueError):
        func.calc(x, -1)
    # Test calc_inv
    output = func.calc_inv(y, 2)
    assert output.shape == (3,)
    assert_allclose(output[0], func.value_inv(y))
    assert_allclose(output[1], func.deriv_inv(y))
    assert_allclose(output[2], func.deriv2_inv(y))
    output = func.calc_inv(y, 1)
    assert output.shape == (2,)
    assert_allclose(output[0], func.value_inv(y))
    assert_allclose(output[1], func.deriv_inv(y))
    output = func.calc_inv(y, 0)
    assert output.shape == (1,)
    assert_allclose(output[0], func.value_inv(y))
    with assert_raises(ValueError):
        func.calc_inv(y, 3)
    with assert_raises(ValueError):
        func.calc_inv(y, -1)
    # Test convenience functions with some chain rules
    assert_allclose(func.value(x), y)
    assert_allclose(func.value_inv(y), x)
    assert_allclose(func.deriv(x)*func.deriv_inv(y), 1.0)
    assert_allclose(func.deriv2_inv(y)*func.deriv(x)**2 + func.deriv_inv(y)*func.deriv2(x),
                    0.0, atol=1e-10)
    assert_allclose(func.deriv2(x)*func.deriv_inv(y)**2 + func.deriv(x)*func.deriv2_inv(y),
                    0.0, atol=1e-10)
    # Test broadcasting calc and calc_inv
    for method in func.calc, func.calc_inv:
        for nderiv in 0, 1, 2:
            assert method(x, nderiv).shape == (nderiv + 1,)
            result = method(np.array([x, 2*x]), nderiv)
            assert result.shape == (2, nderiv + 1)
            assert (result == np.array([method(x, nderiv), method(2*x, nderiv)])).all()
    # Test broadcasting convenience functions
    methods = func.value, func.value_inv, func.deriv, func.deriv_inv, func.deriv2, func.deriv2_inv
    for method in methods:
        assert isinstance(method(x), float)
        result = method(np.array([x, 2*x]))
        assert result.shape == (2,)
        assert (result == np.array([method(x), method(2*x)])).all()
