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
from qcgrids.ext import Cellgrid, Exp, Ln, Linear, Identity, Constant, Power, Rational, \
    UniformCubicSpline, Composed


def test_something():
    cell = Cell(np.identity(3)*5.0)
    grid = Cellgrid(cell)
    grid.append_many(np.random.uniform(0, 10, (2, 3)), np.random.uniform(0, 1, 2))
    assert grid.npoint == 2


def test_exp():
    func = Exp(0.2, -0.3)
    assert func.slope == 0.2
    assert func.offset == -0.3
    check_fn(func, True, 1.2, np.exp(0.2*1.2 - 0.3))


def test_ln():
    func = Ln(0.2, 0.3)
    assert func.prefac == 0.2
    assert func.alpha == 0.3
    check_fn(func, True, 1.2, 0.2*np.log(0.3*1.2))


def test_linear():
    func = Linear(0.2, 0.3)
    assert func.slope == 0.2
    assert func.offset == 0.3
    check_fn(func, True, 1.2, 0.2*1.2 + 0.3)


def test_identity():
    func = Identity()
    check_fn(func, True, 1.2, 1.2)


def test_constant():
    func = Constant(0.3)
    check_fn(func, False, 1.2, 0.3)


def test_power():
    func = Power(0.3, 2.1)
    assert func.prefac == 0.3
    assert func.power == 2.1
    check_fn(func, True, 1.2, 0.3*1.2**2.1)


def test_rational():
    func = Rational(0.3, 2.1)
    assert func.prefac == 0.3
    assert func.root == 2.1
    check_fn(func, True, 1.2, 0.3*1.2/(1.0 - 1.2/2.1))


def test_uniform_cubic_spline1():
    values = np.array([1.2, 0.3, -4.0])
    derivs = np.array([-0.3, -1.0, 0.2])
    func = UniformCubicSpline(values, derivs)
    assert func.npoint == 3
    assert (func.values == values).all()
    assert (func.derivs == derivs).all()
    assert func.derivs2 is None
    assert func.left == 0.0
    assert func.right == 2.0
    assert func.x(0) == 0.0
    assert func.x(1) == 1.0
    assert func.x(2) == 2.0
    check_fn(func, False, 1.2, -0.2816)


def test_uniform_cubic_spline2():
    values = np.array([1.2, 0.3, -4.0])
    func = UniformCubicSpline(values)
    assert func.npoint == 3
    assert (func.values == values).all()
    assert (func.derivs != 0).any()
    assert func.derivs2 is None
    assert func.left == 0.0
    assert func.right == 2.0
    assert func.x(0) == 0.0
    assert func.x(1) == 1.0
    assert func.x(2) == 2.0
    check_fn(func, False, 1.2, -0.3152)
    derivs = func.derivs
    derivs_backup = func.derivs.copy()
    func.derivs[:] = 0.0
    assert (derivs == 0.0).all()
    func.fit_derivs()
    assert_allclose(derivs, derivs_backup)
    assert_allclose(func.derivs, derivs_backup)


def test_composed1():
    exp1 = Exp(0.5, 1.5)
    exp2 = Exp(0.1, 2.5)
    linear = Linear(0.4, 1.2)
    identity = Identity()
    values = np.array([0.5, 1.0, 2.3])
    derivs = np.array([-0.1, 1.2, 2.0])
    composed = Composed(exp1, exp2, linear, identity, values, derivs)
    assert composed.spline.npoint == 3
    assert composed.x_transform is exp1
    assert composed.y_transform is exp2
    assert composed.left_extra is linear
    assert composed.right_extra is identity
    check_fn(composed, False, 10.0, 1.2192817112639593)
    check_fn(composed, False, -0.5, 0.4*(-0.5) + 1.2)
    check_fn(composed, False, 15.5, 15.5)


def test_composed2():
    exp1 = Exp(0.5, 1.5)
    exp2 = Exp(0.1, 2.5)
    linear = Linear(0.4, 1.2)
    identity = Identity()
    values = np.array([0.5, 1.0, 2.3])
    composed = Composed(exp1, exp2, linear, identity, values)
    check_fn(composed, False, 10.0, 1.6362496412909484)
    check_fn(composed, False, -0.5, 0.4*(-0.5) + 1.2)
    check_fn(composed, False, 15.5, 15.5)


# pylint: disable=too-many-statements
def check_fn(func, invertible, x, y):
    """Test if a subclass of ScalarFunction has the right API from the base class.

    Parameters
    ----------
    func : ScalarFunction
        Function to be tested.
    invertible : bool
        The expected value of the invertible attribute of func.
    x : float
        The x-value to test at.
    y : float
        The y-value to test at.

    """
    # Test invertible
    assert func.invertible == invertible
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
    if invertible:
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
    if invertible:
        assert_allclose(func.value_inv(y), x)
        assert_allclose(func.deriv(x)*func.deriv_inv(y), 1.0)
        assert_allclose(func.deriv2_inv(y)*func.deriv(x)**2 + func.deriv_inv(y)*func.deriv2(x),
                        0.0, atol=1e-10)
        assert_allclose(func.deriv2(x)*func.deriv_inv(y)**2 + func.deriv(x)*func.deriv2_inv(y),
                        0.0, atol=1e-10)
    # Test broadcasting calc and calc_inv
    if invertible:
        methods = func.calc, func.calc_inv
    else:
        methods = func.calc,
    for method in methods:
        for nderiv in 0, 1, 2:
            assert method(x, nderiv).shape == (nderiv + 1,)
            result = method(np.array([x, 2*x]), nderiv)
            assert result.shape == (2, nderiv + 1)
            assert issubclass(result.dtype.type, float)
            assert (result == np.array([method(x, nderiv), method(2*x, nderiv)])).all()
            result = method(np.array([int(x), 2*int(x) + 1]), nderiv)
            assert result.shape == (2, nderiv + 1)
            assert issubclass(result.dtype.type, float)
            assert (result == np.array([
                method(int(x), nderiv), method(2*int(x) + 1, nderiv)])).all()
    # Test broadcasting convenience functions
    if invertible:
        methods = (func.value, func.value_inv, func.deriv, func.deriv_inv, func.deriv2,
                   func.deriv2_inv)
    else:
        methods = func.value, func.deriv, func.deriv2
    for method in methods:
        assert isinstance(method(x), float)
        result = method(np.array([x, 2*x]))
        assert result.shape == (2,)
        assert (result == np.array([method(x), method(2*x)])).all()
        result = method(np.array([int(x), 2*int(x) + 1]))
        assert result.shape == (2,)
        assert issubclass(result.dtype.type, float)
        assert (result == np.array([method(int(x)), method(2*int(x) + 1)])).all()
