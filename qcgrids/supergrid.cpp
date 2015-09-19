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
#include <memory>
#include <vector>

#include <celllists/cell.h>
#include <celllists/decomposition.h>
#include <celllists/iterators.h>
#include <celllists/vec3.h>

#include "qcgrids/supergrid.h"


namespace cl = celllists;
namespace vec3 = celllists::vec3;

namespace qcgrids {


SupergridPoint::SupergridPoint(const double* cart, const double weight,
    const size_t index) : icell_{0, 0, 0}, weight_(weight), index_(index) {
  std::copy(cart, cart + 3, cart_);
}


Supergrid::Supergrid(const cl::Cell& cell, double spacing)
    : grid_array_(),
      cell_(new cl::Cell(cell)),
      subcell_(nullptr),
      shape_{0, 0, 0},
      cell_map_(nullptr) {
  subcell_.reset(cell_->create_subcell(spacing, shape_));
}


const cl::CellMap* Supergrid::cell_map() const {
  if (cell_map_.get() == nullptr)
    throw std::logic_error("cell_map can not be requested before calling sort.");
  return cell_map_.get();
}


void Supergrid::emplace_back(const double* cart, const double weight) {
  grid_array_.emplace_back(cart, weight, grid_array_.size());
}


void Supergrid::emplace_back_many(const double* cart, const double* weight,
  const size_t npoint) {
  for (size_t ipoint = 0; ipoint < npoint; ++ipoint) {
    grid_array_.emplace_back(cart, *weight, grid_array_.size());
    cart += 3;
    ++weight;
  }
}


void Supergrid::sort() {
  // Do domain decomposition
  cl::assign_icell(*subcell_, shape_, grid_array_.data(), grid_array_.size(),
      sizeof(SupergridPoint));
  cl::sort_by_icell(grid_array_.data(), grid_array_.size(), sizeof(SupergridPoint));
  cell_map_.reset(cl::create_cell_map(grid_array_.data(), grid_array_.size(),
      sizeof(SupergridPoint)));
}



Subgrid* Supergrid::create_subgrid(const double* center, const double cutoff) const {
  // The sort method must have been called before
  if (cell_map_.get() == nullptr)
    throw std::logic_error("sort must be called before calling create_subgrid.");

  // Use domains to to this efficiently
  std::vector<int> bars;
  subcell_->bars_cutoff(center, cutoff, &bars);
  Subgrid* subgrid = new Subgrid(center);
  for (cl::BarIterator3D bit(bars, shape_); bit.busy(); ++bit) {
    auto it = cell_map_->find(bit.icell());
    if (it != cell_map_->end()) {
      for (size_t ipoint = it->second[0]; ipoint < it->second[1]; ++ipoint) {
        double cart[3];
        std::copy(grid_array_[ipoint].cart_, grid_array_[ipoint].cart_ + 3, cart);
        cell_->iadd_vec(cart, bit.coeffs());
        double d = vec3::distance(center, cart);
        if (d < cutoff)
          subgrid->emplace_back(cart, d, grid_array_[ipoint].weight_, ipoint);
      }
    }
  }

  // Sort by distance from the center and return.
  subgrid->sort();
  return subgrid;
}



}  // namespace qcgrids

// vim: textwidth=90 et ts=2 sw=2
