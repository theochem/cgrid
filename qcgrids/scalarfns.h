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


#include <memory>


#ifndef QCGRIDS_SCALARFNS_H_
#define QCGRIDS_SCALARFNS_H_


namespace qcgrids {

/** @brief
        Abstract base class for (parameterized) scalar functions.

    In the definition of integration grids, we often need to specify transformations and
    other functions, which may contain parameters. These are objects of the ScalarFunction
    class, or its derivatives.

    Derived classes must override ::calc and ::calc_inv. They may also override all other
    methods if they wish to make these functions more efficient. For simple functions, the
    latter is recommended.
*/
class ScalarFunction {
 public:
  const bool invertible;  //!< Flag for invertible functions.

  //! Create a ScalarFunction object
  explicit ScalarFunction(bool invertible) : invertible(invertible) {}
  virtual ~ScalarFunction() {}

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


  // Convenience functions

  /** @brief
        Compute the value of the function

      param x
        The input to the function.

      @return
        The function value.
  */
  virtual double value(const double x) const;

  /** @brief
        Compute the value of the inverse function

      param y
        The input to the inverse function.

      @return
        The inverse function value.
  */
  virtual double value_inv(const double y) const;

  /** @brief
        Compute the derivative of the function

      param x
        The input to the function.

      @return
        The function value.
  */
  virtual double deriv(const double x) const;

  /** @brief
        Compute the derivative of the inverse function

      param y
        The input to the inverse function.

      @return
        The inverse function value.
  */
  virtual double deriv_inv(const double y) const;

  /** @brief
        Compute the second derivative of the function

      param x
        The input to the function.

      @return
        The function value.
  */
  virtual double deriv2(const double x) const;

  /** @brief
        Compute the second derivative of the inverse function

      param y
        The input to the inverse function.

      @return
        The inverse function value.
  */
  virtual double deriv2_inv(const double y) const;
};


//! An exponential function.
class Exp : public ScalarFunction {
 public:
  const double prefac;  //!< prefac in y=prefac*exp(alpha*x)
  const double alpha;   //!< alpha in y=prefac*exp(alpha*x)

  Exp() = delete;
  //! Create an instance of Exp.
  explicit Exp(double prefac, double alpha) : ScalarFunction(true), prefac(prefac),
      alpha(alpha) {}
  virtual ~Exp() {}

  void calc(const double x, const int nderiv, double* const output) const;
  void calc_inv(const double y, const int nderiv, double* const output) const;
  virtual double value(const double x) const;
  virtual double value_inv(const double y) const;
  virtual double deriv(const double x) const;
  virtual double deriv_inv(const double x) const;
  virtual double deriv2(const double x) const;
  virtual double deriv2_inv(const double x) const;
};


//! A logarithmic function.
class Ln : public ScalarFunction {
 public:
  const double prefac;  //!< The prefac, A in y=A*ln(alpha*x)
  const double alpha;   //!< The alpha, alpha in y=A*ln(alpha*x)

  Ln() = delete;
  //! Create an instance of Exp with given alpha.
  explicit Ln(double prefac, double alpha) : ScalarFunction(true), prefac(prefac),
    alpha(alpha) {}
  virtual ~Ln() {}

  void calc(const double x, const int nderiv, double* const output) const;
  void calc_inv(const double y, const int nderiv, double* const output) const;
  virtual double value(const double x) const;
  virtual double value_inv(const double y) const;
  virtual double deriv(const double x) const;
  virtual double deriv_inv(const double x) const;
  virtual double deriv2(const double x) const;
  virtual double deriv2_inv(const double x) const;
};



/** @brief
      An abstract base class for all splines.
*/
class Spline : public ScalarFunction {
 public:
  const size_t npoint;  //!< The number of points in the spline.

  Spline() = delete;
  //! Create a splin with npoint points.
  explicit Spline(size_t npoint) : ScalarFunction(false), npoint(npoint) {}

  //! The position of the left-most grid point.
  virtual double left() const = 0;
  //! The position of the right-most grid point.
  virtual double right() const = 0;
  //! Return position on x-axis of point with given index.
  virtual double x(const size_t index) const = 0;

  void calc_inv(const double y, const int nderiv, double* const output) const;
};


/** @brief
      Simple solver for a symmetric tridiagonal system of linear equations.

    rhs and diag2 are modified in-pace.

    @param diag1
      The diagonal of the matrix. (n elements)

    @param diag2
      The second diagonal. (n-1 elements)

    @param rhs
      The right-hand side. (n elements)

    @param solution
      The output vector. (n elements)

    @param n
      The number of equations and unknowns.
*/
void tridiagsym_solve(double* diag1, double* diag2, double* rhs, double* solution, size_t n);


/** @brief
      A cubic spline implementation, with fixed uniformly distributed grid points 0, 1, 2,
      ..., npoint-1.
*/
class UniformCubicSpline : public Spline {
 public:
  // Not using std::vector to keep things simple: constant pointer to non-constant data.
  double * const values;  //!< Array with the function value at the grid points.
  double * const derivs;  //!< Array with the derivative at the grid points.

  //! Construct new spline with npoints, zeroing values and derivs.
  explicit UniformCubicSpline(size_t npoint) : Spline(npoint), values(new double[npoint]),
      derivs(new double[npoint]) {}
  //! Create new spline using function values at grid point. Derivatives are fitted.
  explicit UniformCubicSpline(size_t npoint, const double* const values);
  //! Create new spline using function values and derivatives at grid point.
  explicit UniformCubicSpline(size_t npoint, const double* const values,
      const double* const derivs);
  //! Copy constructor
  UniformCubicSpline(const UniformCubicSpline &obj);
  //! Copy assignment is not always possible because npoint is constant.
  UniformCubicSpline& operator=(const UniformCubicSpline &obj) = delete;

  virtual ~UniformCubicSpline();

  double left() const { return 0.0; }
  double right() const { return static_cast<double>(npoint) - 1.0; }
  double x(const size_t index) const { return static_cast<double>(index); }

  void calc(const double x, const int nderiv, double* const output) const;

  /** @brief
        Determine the derivatives at the grid points for the case of a natural spline.

      A natural cubic splint has vanishing second derivatives are edges and continuous
      second derivatives at all other points.
  */
  void fit_derivs();
};


}  // namespace qcgrids

#endif  // QCGRIDS_SCALARFNS_H_

// vim: textwidth=90 et ts=2 sw=2
