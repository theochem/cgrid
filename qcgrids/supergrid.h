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

/** @file */


#ifndef QCGRIDS_SUPERGRID_H_
#define QCGRIDS_SUPERGRID_H_


#include <memory>
#include <vector>

#include <celllists/cell.h>
#include <celllists/decomposition.h>

#include "qcgrids/subgrid.h"


namespace cl = celllists;

namespace qcgrids {


class SupergridPoint {
 public:
  SupergridPoint(const double* cart, const double weight, const size_t index);
  SupergridPoint() : cart_{0.0, 0.0, 0.0}, icell_{0, 0, 0}, weight_(0.0), index_(0) {}

  double cart_[3];
  int icell_[3];
  double weight_;
  size_t index_;
};


class Supergrid {
 public:
  explicit Supergrid(const cl::Cell& cell, double spacing = 1.0);

  const std::vector<SupergridPoint>& grid_array() const { return grid_array_; }
  const cl::Cell* cell() const { return cell_.get(); }
  const cl::Cell* subcell() const { return subcell_.get(); }
  const int* shape() const { return shape_; }
  const cl::CellMap* cell_map() const;

  void emplace_back(const double* cart, const double weight);
  void emplace_back_many(const double* cart, const double* weight, const size_t npoint);
  void sort();
  Subgrid* create_subgrid(const double* center, const double cutoff) const;

 private:
  std::vector<SupergridPoint> grid_array_;
  std::unique_ptr<cl::Cell> cell_;
  std::unique_ptr<cl::Cell> subcell_;
  int shape_[3];
  std::unique_ptr<cl::CellMap> cell_map_;
};


}  // namespace qcgrids


#endif  // QCGRIDS_SUPERGRID_H_

// vim: textwidth=90 et ts=2 sw=2
