// QCGrids is a numerical integration library for quantum chemistry.
// Copyright (C) 2011-2015 Toon Verstraelen
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
//--


#include "qcgrids/subgrid.h"

#include <algorithm>
#include <stdexcept>


namespace qcgrids {


SubgridPoint::SubgridPoint(const double* cart, const double distance, const double weight,
    const int index) : distance_(distance), weight_(weight), index_(index) {
  std::copy(cart, cart + 3, cart_);
}


Subgrid::Subgrid(const double* center) : grid_array_() {
  std::copy(center, center + 3, center_);
}


void Subgrid::iadd_super(double* super_array, const double* sub_array) {
  size_t isub = 0;
  for (const SubgridPoint& point : grid_array_) {
    super_array[point.index_] += sub_array[isub];
    ++isub;
  }
}


void Subgrid::take_sub(const double* super_array, double* sub_array) {
  size_t isub = 0;
  for (const SubgridPoint& point : grid_array_) {
    sub_array[isub] = super_array[point.index_];
    ++isub;
  }
}


double Subgrid::integrate(const double*const* factors, const size_t* strides,
    const size_t nfactor) {
  size_t isub = 0;
  double result = 0.0;
  for (const SubgridPoint& point : grid_array_) {
    double term = point.weight_;
    for (size_t ifactor = 0; ifactor < nfactor; ++ifactor) {
      term *= factors[ifactor][isub*strides[ifactor]];
    }
    result += term;
    ++isub;
  }
  return result;
}


void Subgrid::emplace_back(const double* cart, const double distance, const double weight,
    const int index) {
  grid_array_.emplace_back(cart, distance, weight, index);
}


void Subgrid::sort() {
  std::sort(grid_array_.begin(), grid_array_.end());
}



}  // namespace qcgrids

// vim: textwidth=90 et ts=2 sw=2
