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



#include "qcgrids/scalarfns.h"

#include <stdexcept>

#include <cmath>


namespace qcgrids {


double ScalarFunction::calc(const double x) const {
  double result;
  calc(x, 0, &result);
  return result;
}


double ScalarFunction::calc_inv(const double y) const {
  double result;
  calc_inv(y, 0, &result);
  return result;
}


void Exponential::calc(const double x, const int nderiv, double* const output) const {
  if (nderiv < 0) throw std::domain_error("nderiv cannot be negative.");
  output[0] = exp(exponent*x);
  if (nderiv > 0) output[1] = output[0]*exponent;
  if (nderiv > 1) output[2] = output[1]*exponent;
  if (nderiv > 2) throw std::domain_error("nderiv cannot be larger than two.");
}


void Exponential::calc_inv(const double y, const int nderiv, double* const output) const {
  if (nderiv < 0) throw std::domain_error("nderiv cannot be negative.");
  output[0] = log(y)/exponent;
  if (nderiv > 0) output[1] = 1.0/(y*exponent);
  if (nderiv > 1) output[2] = -output[1]/y;
  if (nderiv > 2) throw std::domain_error("nderiv cannot be larger than two.");
}


}  // namespace qcgrids

// vim: textwidth=90 et ts=2 sw=2
