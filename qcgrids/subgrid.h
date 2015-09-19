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


#ifndef QCGRIDS_SUBGRID_H_
#define QCGRIDS_SUBGRID_H_

#include <cstdlib>
#include <vector>


namespace qcgrids {


class SubgridPoint {
 public:
  SubgridPoint(const double* cart, const double distance, const double weight,
      const size_t index);
  SubgridPoint() : cart_{0.0, 0.0, 0.0}, distance_(0.0), weight_(0), index_(0) {}

  bool operator<(const SubgridPoint& other) const { return distance_ < other.distance_; }

  double cart_[3];
  double distance_;
  double weight_;
  size_t index_;
};


class Subgrid {
 public:
  explicit Subgrid(const double* center);

  std::vector<SubgridPoint>& grid_array() { return grid_array_; }
  const double* center() { return center_; }

  void iadd_super(double* super_array, const double* sub_array);
  void take_sub(const double* super_array, double* sub_array);
  double integrate(const double*const* factors, const size_t* strides,
      const size_t nfactor);

  void emplace_back(const double* cart, const double distance, const double weight,
      const size_t index);
  void sort();

 private:
  std::vector<SubgridPoint> grid_array_;
  double center_[3];
};


}  // namespace qcgrids


#endif  // QCGRIDS_SUBGRID_H_

// vim: textwidth=90 et ts=2 sw=2
