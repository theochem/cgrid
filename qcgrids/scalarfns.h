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

/** @file
    @brief
      Parameterized scalar functions, their inverses and first and second derivatives.

    These scalar functions, including splines, are an important auxiliary for the
    definition of several different types of numerical integration grids. A special case
    is the composed function, which can be used to represent radial density or potential
    functions with splines, including the correct functional behavior towards the limits
    (defined by extrapolation functions).
*/


#include <memory>

#include <cmath>


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
  //! Create a ScalarFunction object
  explicit ScalarFunction(bool invertible) : invertible_(invertible) {}
  virtual ~ScalarFunction() {}

  //! Return the invertible flag. Only true when function is guaranteed to be invertible.
  bool invertible() const { return invertible_; }

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
  virtual void calc_inv(const double y, const int nderiv, double* const output) const;


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

 private:
  const bool invertible_;  //!< Flag for invertible functions.
};


//! An exponential function.
class Exp : public ScalarFunction {
 public:
  Exp() = delete;
  //! Create an instance of Exp.
  explicit Exp(double slope, double offset);
  virtual ~Exp() {}

  double slope() const { return slope_; }    //!< slope in y = exp(slope*x + offset)
  double offset() const { return offset_; }  //!< offset in y = exp(slope*x + offset)

  void calc(const double x, const int nderiv, double* const output) const;
  void calc_inv(const double y, const int nderiv, double* const output) const;
  double value(const double x) const { return exp(slope_*x + offset_); }
  double value_inv(const double y) const { return (log(y) - offset_)/slope_; }
  double deriv(const double x) const { return slope_*exp(slope_*x + offset_); }
  double deriv_inv(const double y) const { return 1.0/(y*slope_); }
  double deriv2(const double x) const { return slope_*slope_*exp(slope_*x + offset_); }
  double deriv2_inv(const double y) const { return -1.0/(y*y*slope_); }

 private:
  const double slope_;   //!< slope in y = exp(slope*x + offset)
  const double offset_;  //!< offset in y = exp(slope*x + offset)
};


//! A logarithmic function.
class Ln : public ScalarFunction {
 public:
  Ln() = delete;
  //! Create an instance of Ln.
  explicit Ln(double prefac, double alpha);
  virtual ~Ln() {}

  double prefac() const { return prefac_; }  //!< prefac in y = prefac*ln(alpha*x)
  double alpha() const { return alpha_; }    //!< alpha in y = prefac*ln(alpha*x)

  void calc(const double x, const int nderiv, double* const output) const;
  void calc_inv(const double y, const int nderiv, double* const output) const;
  double value(const double x) const { return prefac_*log(alpha_*x); }
  double value_inv(const double y) const { return exp(y/prefac_)/alpha_; }
  double deriv(const double x) const { return prefac_/x; }
  double deriv_inv(const double y) const { return exp(y/prefac_)/(alpha_*prefac_); }
  double deriv2(const double x) const { return -prefac_/(x*x); }
  double deriv2_inv(const double y) const { return exp(y/prefac_)/(alpha_*prefac_*prefac_); }

 private:
  const double prefac_;  //!< prefac in y = prefac*ln(alpha*x)
  const double alpha_;   //!< alpha in y = prefac*ln(alpha*x)
};


//! A linear function.
class Linear : public ScalarFunction {
 public:
  Linear() = delete;
  //! Create an instance of Ln.
  explicit Linear(double slope, double offset);
  virtual ~Linear() {}

  double slope() const { return slope_; }     //!< slope in y = slope*x + offset
  double offset() const { return offset_; }   //!< offset in y = slope*x + offset

  void calc(const double x, const int nderiv, double* const output) const;
  void calc_inv(const double y, const int nderiv, double* const output) const;
  double value(const double x) const { return slope_*x + offset_; }
  double value_inv(const double y) const { return (y - offset_)/slope_; }
  double deriv(const double) const { return slope_; }
  double deriv_inv(const double) const { return 1.0/slope_; }
  double deriv2(const double) const { return 0.0; }
  double deriv2_inv(const double) const { return 0.0; }

 private:
  const double slope_;    //!< slope in y = slope*x + offset
  const double offset_;   //!< offset in y = slope*x + offset
};


//! Identity function.
class Identity : public ScalarFunction {
 public:
  //! Create an instance of Identity.
  Identity() : ScalarFunction(true) {}
  virtual ~Identity() {}

  void calc(const double x, const int nderiv, double* const output) const;
  void calc_inv(const double y, const int nderiv, double* const output) const;
  double value(const double x) const { return x; }
  double value_inv(const double y) const { return y; }
  double deriv(const double) const { return 1.0; }
  double deriv_inv(const double) const { return 1.0; }
  double deriv2(const double) const { return 0.0; }
  double deriv2_inv(const double) const { return 0.0; }
};


//! Constant function.
class Constant : public ScalarFunction {
 public:
  Constant() = delete;
  //! Create an instance of Identity.
  explicit Constant(double offset) : ScalarFunction(false), offset_(offset) {}
  virtual ~Constant() {}

  double offset() const { return offset_; }  //!< Value of the constant function.

  void calc(const double x, const int nderiv, double* const output) const;
  double value(const double) const { return offset_; }
  double deriv(const double) const { return 0.0; }
  double deriv2(const double) const { return 0.0; }

 private:
  const double offset_;  //!< Value of the constant function.
};


//! A Power function, only for x > 0.
class Power : public ScalarFunction {
 public:
  Power() = delete;
  /** @brief
        Create an instance of Power.

      Note: The function is marked as invertible but this only true when the function is
      restricted to positive x. This is what we need in practice, so we'll imply this
      limitation. (No strict checking.)
  */
  explicit Power(double prefac, double power);
  virtual ~Power() {}

  double prefac() const { return prefac_; }  //!< prefac in y = prefac*(x**power)
  double power() const { return power_; }    //!< power in y = prefac*(x**power)

  void calc(const double x, const int nderiv, double* const output) const;
  void calc_inv(const double y, const int nderiv, double* const output) const;
  double value(const double x) const { return prefac_*pow(x, power_); }
  double value_inv(const double y) const { return pow(y/prefac_, 1.0/power_); }
  double deriv(const double x) const { return prefac_*power_*pow(x, power_ - 1); }
  double deriv_inv(const double y) const
    { return pow(y/prefac_, 1.0/power_ - 1.0)/(power_*prefac_); }
  double deriv2(const double x) const;
  double deriv2_inv(const double y) const;

 private:
  const double prefac_;  //!< prefac in y = prefac*(x**power)
  const double power_;   //!< power in y = prefac*(x**power)
};


//! A basic rational function: y(x) = prefac*x/(1 - x/root), not for x == root.
class Rational : public ScalarFunction {
 public:
  Rational() = delete;
  //! Create an instance of Rational.
  explicit Rational(double prefac, double root);
  virtual ~Rational() {}

  double prefac() const { return prefac_; }  //!< prefac in prefac*x/(1 - x/root)
  double root() const { return root_; }      //!< root in prefac*x/(1 - x/root)

  void calc(const double x, const int nderiv, double* const output) const;
  void calc_inv(const double y, const int nderiv, double* const output) const;
  double value(const double x) const { return prefac_*x/(1 - x/root_); }
  double value_inv(const double y) const { return y/(prefac_ + y/root_); }
  double deriv(const double x) const
    { double t = 1 - x/root_; return prefac_/(t*t);  }
  double deriv_inv(const double y) const
    { double t = prefac_ + y/root_; return prefac_/(t*t); }
  double deriv2(const double x) const
    { double t = 1 - x/root_; return 2*prefac_/(root_*t*t*t); }
  double deriv2_inv(const double y) const
    { double t = prefac_ + y/root_; return -2*prefac_/(root_*t*t*t); }

 private:
  const double prefac_;  //!< prefac in prefac*x/(1 - x/root)
  const double root_;    //!< root in prefac*x/(1 - x/root)
};



/** @brief
      An abstract base class for all splines.
*/
class Spline : public ScalarFunction {
 public:
  Spline() = delete;
  //! Create a spline with npoint points.
  explicit Spline(size_t npoint);

  size_t npoint() const { return npoint_; }  //!< The number of points in the spline.

  // Derived classes may return NULL for any of the following three functions.
  virtual double* values() const = 0;   //!< Array with the function value at the grid points.
  virtual double* derivs() const = 0;   //!< Array with the derivative at the grid points.
  virtual double* derivs2() const = 0;  //!< Array with the second derivative at the grid points.

  //! The position of the left-most grid point.
  virtual double left() const = 0;
  //! The position of the right-most grid point.
  virtual double right() const = 0;
  //! Return position on x-axis of point with given index.
  virtual double x(const size_t index) const = 0;

  /** @brief
        Determine the (second) derivatives at the grid points for the case of a natural
        spline.

      A natural cubic splint has vanishing second derivatives are edges and continuous
      second derivatives at all other points.
  */
  virtual void fit_derivs() = 0;

 protected:
  const size_t npoint_;  //!< The number of points in the spline.
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
  //! Construct new spline with npoints, zeroing values and derivs.
  explicit UniformCubicSpline(size_t npoint);
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

  double* values() const { return values_; }
  double* derivs() const { return derivs_; }
  double* derivs2() const { return NULL; }

  double left() const { return 0.0; }
  double right() const { return static_cast<double>(npoint_) - 1.0; }
  double x(const size_t index) const { return static_cast<double>(index); }

  void calc(const double x, const int nderiv, double* const output) const;

  void fit_derivs();

 private:
  // Not using std::vector to keep things simple: constant pointer to non-constant data.
  double * const values_;  //!< Array with the function value at the grid points.
  double * const derivs_;  //!< Array with the derivative at the grid points.
};


/** @brief
      Combination of spline and different functions for extrapolation and transformation.

    Two functions can be provided for extrapolation on the left and the right of the
    spline. Also a transformation of the input and output can be given. Transformation
    functions must be invertible.

    The composed function has the form y_transf(spline(inv_x_transf(x))) when
    inv_x_transf(x) falls in the interval for which the underlying spline is define. The
    extrapolation functions are called otherwise.
*/
class Composed : public ScalarFunction {
 public:
  Composed() = delete;
  /**  @brief
        Create a composed function based on a spline and some other functions. Spline will
        not be modified.
  */
  explicit Composed(Spline* spline, ScalarFunction* x_transform, ScalarFunction* y_transform,
                    ScalarFunction* left_extra, ScalarFunction* right_extra);

  /**  @brief
        Create a composed function based on a spline and some other functions. Spline will
        derivatives will be fitted to inverse y-transform of values.
  */
  explicit Composed(Spline* spline, ScalarFunction* x_transform, ScalarFunction* y_transform,
                    ScalarFunction* left_extra, ScalarFunction* right_extra,
                    const double* values);

  /**  @brief
        Create a composed function based on a spline and some other functions. Spline will
        be configured by inverse y-transformation of values and derivs.
  */
  explicit Composed(Spline* spline, ScalarFunction* x_transform, ScalarFunction* y_transform,
                    ScalarFunction* left_extra, ScalarFunction* right_extra,
                    const double* values, const double* derivs);

  virtual ~Composed();

  Spline* spline() const { return spline_; }                    //!< The underlying cubic spline.
  ScalarFunction* x_transform() const { return x_transform_; }  //!< The transformation for x.
  ScalarFunction* y_transform() const { return y_transform_; }  //!< The transformation for y.
  ScalarFunction* left_extra() const { return left_extra_; }    //!< The left extrapolation.
  ScalarFunction* right_extra() const { return right_extra_; }  //!< The right extrapolation.

  void calc(const double x, const int nderiv, double* const output) const;

 private:
  Identity* const default_transform_;  //!< The default transformation for x and y.
  Constant* const default_extra_;      //!< The default extrapolation for left and right.
  Spline* const spline_;               //!< The underlying cubic spline.
  ScalarFunction* const x_transform_;  //!< The transformation for x.
  ScalarFunction* const y_transform_;  //!< The transformation for y.
  ScalarFunction* const left_extra_;   //!< The left extrapolation.
  ScalarFunction* const right_extra_;  //!< The right extrapolation.
};


}  // namespace qcgrids


#endif  // QCGRIDS_SCALARFNS_H_


// vim: textwidth=90 et ts=2 sw=2
