// CGrid is a library for for molecular numerical integration.
// Copyright (C) 2011-2019 The CGrid Development Team
//
// This file is part of CGrid.
//
// CGrid is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 3
// of the License, or (at your option) any later version.
//
// CGrid is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, see <http://www.gnu.org/licenses/>
// --


#ifndef CGRID_TESTS_RIDDERS_H_
#define CGRID_TESTS_RIDDERS_H_


#include <algorithm>
#include <stdexcept>
#include <vector>

/** @brief
      Approximate the derivative of a scalar function with Ridders' method.

    @param function
      The function to be differentiated.

    @param origin
      The point where the derivative must be computed.

    @param stepsize
      The (initial) step size for the finite difference calculation.

    @param con
      The rate at which the step size is decreased (contracted). Must be larger than one.

    @param safe
      The safety check used to terminate the algorithm. If Errors between successive
      orders become larger than ``safe`` times the error on the best estimate, the
      algorithm stop. This happens due to round-off errors.

    @param maxiter
      The maximum number of iterations, equals the maximum number of function calls and
      also the highest polynomial order in the Neville method.

    @param outerror
      Optional output argument with the error estimate for the derivative.

    @returns
      The derivative.
*/


/*
  TODO: write a few overloads that also work with simpler cases: plain function (instead
  of method) and plain function with extra arguments.
*/

template<typename T, typename FLT>
FLT diff_ridders(
    const T* obj, FLT (T::*function)(const FLT) const, FLT origin, FLT stepsize,
    FLT con = 1.4, FLT safe = 2.0, size_t maxiter = 15, FLT* outerror = NULL) {

  if (stepsize == 0.0)
      throw std::domain_error("stepsize must be nonzero.");
  if (con <= 1.0)
      throw std::domain_error("con must be larger than one.");
  if (safe <= 1.0)
      throw std::domain_error("safe must be larger than one.");

  // Initialization for Neville loop.
  FLT con2 = con*con;
  std::vector<std::vector<FLT>> table = {{
      ((obj->*function)(origin + stepsize) - (obj->*function)(origin - stepsize)) /
      (2.0*stepsize) }};
  FLT estimate = NAN;
  FLT error = NAN;

  // Loop based on Neville's method.
  // Successive rows in the table will go to smaller stepsizes.
  // Successive columns in the table go to higher orders of extrapolation.
  for (size_t i = 0; i < maxiter; i++) {
    // Reduce step size.
    stepsize /= con;
    // First-order approximation at current step size.
    table.push_back({
        ((obj->*function)(origin + stepsize) - (obj->*function)(origin - stepsize)) /
        (2.0*stepsize) });
    // Compute higher-orders
    FLT fac = con2;
    for (size_t j = 1; j <= i; j++) {
      // Compute extrapolations of various orders, requiring no new
      // function evaluations. This is a recursion relation based on
      // Neville's method.
      table.at(i).push_back((table.at(i).at(j-1)*fac - table.at(i-1).at(j-1))/(fac-1.0));
      fac = con2*fac;

      // The error strategy is to compare each new extrapolation to one
      // order lower, both at the present stepsize and the previous one:
      FLT current_error = std::max(fabs(table.at(i).at(j) - table.at(i).at(j-1)),
                                   fabs(table.at(i).at(j) - table.at(i-1).at(j-1)));

      // If error has decreased, save the improved estimate.
      if (std::isnan(error) || (current_error <= error)) {
        error = current_error;
        estimate = table.at(i).at(j);
      }
    }

    // If the highest-order estimate is growing larger than the error on the best
    // estimate, the algorithm becomes numerically instable. Time to quit.
    if ((i > 0) && (fabs(table.at(i).at(i) - table.at(i-1).at(i-1)) >= safe * error))
      break;
  }

  if (outerror != NULL)
    *outerror = error;
  return estimate;
}


#endif  // CGRID_TESTS_RIDDERS_H_


// vim: textwidth=90 et ts=2 sw=2
