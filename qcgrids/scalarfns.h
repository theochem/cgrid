// QCGrids is a numerical integration library for quantum chemistry.
// Copyright (C) 2011-2017 The QCGrids developers
//
// This file is part of QCGrids.
//
// QCGrids is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 3
// of the License, or (at your option) any later version.
//
// QCGrids is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, see <http://www.gnu.org/licenses/>
// --

/** @file */


#ifndef QCGRIDS_SCALARFNS_H_
#define QCGRIDS_SCALARFNS_H_


namespace qcgrids {

/** @brief
        Base class for (parameterized) scalar functions.

    In the definition of integration grids, we often need to specify transformations and
    other functions, which may contain parameters. These are objects of the ScalarFunction
    class, or its derivatives.
*/
class ScalarFunction {
 public:
  /** @brief
        Compute the value of the function and optionally some derivatives.

      @param x
        The input to the function.

      @param nderiv
        Derivatives up to this order will be computed. The lowest allowed value is zero.

      @param output
        Pointer to an output array of sufficient size to store the function value and all
        requested derivatives.
  */
  virtual void calc(const double x, const int nderiv, double* const output) const = 0;

  /** @brief
        Compute the value of the inverse function and optionally some derivatives.

      @param y
        The input to the inverse function.

      @param nderiv
        Derivatives up to this order will be computed. The lowest allowed value is zero.

      @param output
        Pointer to an output array of sufficient size to store the inverse function value
        and all requested derivatives.
  */
  virtual void calc_inv(const double y, const int nderiv, double* const output) const = 0;

  /** @brief
        Compute the value of the function

      @params x
        The input to the function.

      @return
        The function value.
  */
  double calc(const double x) const;

  /** @brief
        Compute the value of the inverse function

      @params y
        The input to the inverse function.

      @return
        The inverse function value.
  */
  double calc_inv(const double y) const;
};


//! An exponential function with a fixed exponent.
class Exponential : public ScalarFunction {
 public:
  //! Create an instance of Exponential with given exponent.
  explicit Exponential(double exponent) : exponent(exponent) {}

  // The following lines are needed to inherit the function overloads.
  // They are hidden by default in C++.
  using ScalarFunction::calc;
  using ScalarFunction::calc_inv;

  void calc(const double x, const int nderiv, double* const output) const;
  void calc_inv(const double y, const int nderiv, double* const output) const;
 private:
  const double exponent;
};


}  // namespace qcgrids

#endif  // QCGRIDS_SCALARFNS_H_

// vim: textwidth=90 et ts=2 sw=2
