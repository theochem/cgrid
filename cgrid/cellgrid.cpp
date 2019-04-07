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


#include "cgrid/cellgrid.h"

#include <algorithm>
#include <memory>
#include <vector>

#include <cellcutoff/cell.h>
#include <cellcutoff/decomposition.h>
#include <cellcutoff/iterators.h>
#include <cellcutoff/vec3.h>


namespace cl = cellcutoff;
namespace vec3 = cellcutoff::vec3;

namespace cgrid {


CellgridPoint::CellgridPoint(const double* cart, const double weight,
    const int index) : icell_{0, 0, 0}, weight_(weight), index_(index) {
  std::copy(cart, cart + 3, cart_);
}


Cellgrid::Cellgrid(const cl::Cell& cell, double spacing)
    : grid_array_(),
      cell_(new cl::Cell(cell)),
      subcell_(nullptr),
      shape_{0, 0, 0},
      cell_map_(nullptr) {
  subcell_.reset(cell_->create_subcell(spacing, shape_));
}


const cl::CellMap* Cellgrid::cell_map() const {
  if (cell_map_.get() == nullptr)
    throw std::logic_error("cell_map can not be requested before calling sort.");
  return cell_map_.get();
}


void Cellgrid::emplace_back(const double* cart, const double weight) {
  grid_array_.emplace_back(cart, weight, grid_array_.size());
}


void Cellgrid::emplace_back_many(const double* cart, const double* weight,
  const size_t npoint) {
  for (size_t ipoint = 0; ipoint < npoint; ++ipoint) {
    grid_array_.emplace_back(cart, *weight, grid_array_.size());
    cart += 3;
    ++weight;
  }
}


void Cellgrid::sort() {
  // Do domain decomposition
  cl::assign_icell(*subcell_, shape_, grid_array_.data(), grid_array_.size(),
      sizeof(CellgridPoint));
  cl::sort_by_icell(grid_array_.data(), grid_array_.size(), sizeof(CellgridPoint));
  cell_map_.reset(cl::create_cell_map(grid_array_.data(), grid_array_.size(),
      sizeof(CellgridPoint)));
}


void Cellgrid::iadd_cutoff(const double* center, const double cutoff, GridFunc grid_func,
    const void* extra_arg, double* output) {
  // The sort method must have been called before
  if (cell_map_.get() == nullptr)
    throw std::logic_error("sort must be called before calling create_subgrid.");

  // Loop over all relevant points and compute
  for (cl::DeltaIterator dit(*subcell_, shape_, center, cutoff, grid_array_.data(),
       grid_array_.size(), sizeof(CellgridPoint), *cell_map_); dit.busy(); ++dit) {
    output[dit.ipoint()] += grid_func(dit.delta(), dit.distance(), extra_arg);
  }
}


double Cellgrid::integrate_cutoff(const double* center, const double cutoff,
    GridFunc grid_func, const void* extra_arg, const double* factor) {
  // The sort method must have been called before
  if (cell_map_.get() == nullptr)
    throw std::logic_error("sort must be called before calling create_subgrid.");

  double result = 0.0;
  // Loop over all relevant points and compute
  for (cl::DeltaIterator dit(*subcell_, shape_, center, cutoff, grid_array_.data(),
       grid_array_.size(), sizeof(CellgridPoint), *cell_map_); dit.busy(); ++dit) {
    double term = grid_array_[dit.ipoint()].weight_;
    if (grid_func != nullptr)
      term *= grid_func(dit.delta(), dit.distance(), extra_arg);
    if (factor != nullptr)
      term *= factor[dit.ipoint()];
    result += term;
  }
  return result;
}


}  // namespace cgrid

// vim: textwidth=90 et ts=2 sw=2
