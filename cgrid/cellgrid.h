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

/** @file

    This module implements a grid class that can perform localized numerical integrals and
    evaluation of functions on grids, such that the cost is only proportional to the
    cutoff volume, and not the entire grid size.

 */


#ifndef CGRID_CELLGRID_H_
#define CGRID_CELLGRID_H_


#include <memory>
#include <vector>

#include <cellcutoff/cell.h>
#include <cellcutoff/decomposition.h>


namespace cl = cellcutoff;

namespace cgrid {


//! An integration grid point wit
class CellgridPoint {
 public:
  /** @brief
          Create a cell grid point.

      @param cart
          The three Cartesian coordinates of the grid point.

      @param weight
          The quadrature weight

      @param index
          The index of the grid point in the original order.
   */
  CellgridPoint(const double* cart, const double weight, const int index);

  //! Create a default cell grid point.
  CellgridPoint() : cart_{0.0, 0.0, 0.0}, icell_{0, 0, 0}, weight_(0.0), index_(0) {}

  double cart_[3];   //!< Three Cartesian coordinates of the grid point
  int icell_[3];     //!< Subcell indices in which the grid point sets.
  double weight_;    //!< Quadrature weight
  int index_;        //!< Index of the grid point in the original order.
};


//! Type definition of a Grid function
typedef double (*GridFunc) (const double*, const double, const void*);


//! @brief Numerical integration/evaluation grid with localized features.
class Cellgrid {
 public:
  /** @brief Create a Cellgrid

      @param cell
          The periodic unit cell of the system

      @param spacing
          The spacing between two planes of subcells.
   */
  explicit Cellgrid(const cl::Cell& cell, double spacing = 1.0);

  //! Array of grid points
  const std::vector<CellgridPoint>& grid_array() const { return grid_array_; }

  //! Periodic unit cell
  const cl::Cell* cell() const { return cell_.get(); }
  //! Subcell used for sorting the grid points.
  const cl::Cell* subcell() const { return subcell_.get(); }
  //! Number of repetitions of the subcell along each cell vectors to obtain the periodic box.
  const int* shape() const { return shape_; }
  //! The dictionary of cells with points in each cell.
  const cl::CellMap* cell_map() const;

  //! Add a grid point to the end of the grid.
  void emplace_back(const double* cart, const double weight);
  //! Add many grid points to the end of the grid.
  void emplace_back_many(const double* cart, const double* weight, const size_t npoint);
  //! Sort the grid points (by cell)
  void sort();

  /** @brief
          Evaluate function on grid points within a cutoff sphere and add to total grid.

      @param center
          The center of the cutoff sphere.

      @param cutoff
          The radius of the cutoff sphere.

      @param grid_func
          The function to be evaluated.

      @param extra_arg
          Additional argument to the grid_func.

      @param output
          Array with grid values (total grid).
    */
  void iadd_cutoff(const double* center, const double cutoff, GridFunc grid_func,
      const void* extra_arg, double* output);

  /** @brief
          Compute a localized integral of a grid function times values on the grid.

      @param center
          The center of the cutoff sphere.

      @param cutoff
          The radius of the cutoff sphere.

      @param grid_func
          The function to be evaluated. This contributes a factor to the total integrand.
          When nullptr, it is not evaluated and 1 is used instead.

      @param extra_arg
          Additional argument to the grid_func.

      @param factor
          Array with grid values to be used as the second factor in the integrand. When
          nullptr is giving, this factor is replaced by 1.

      @return
          The value of the numerical integral
    */
  double integrate_cutoff(const double* center, const double cutoff, GridFunc grid_func,
      const void* extra_arg, const double* factor);

 private:
  std::vector<CellgridPoint> grid_array_;   //!< Array with grid points.
  std::unique_ptr<cl::Cell> cell_;          //!< Periodic unit cell.
  std::unique_ptr<cl::Cell> subcell_;       //!< Subcell used for sorting points.
  int shape_[3];                            //!< Repetitions of subcell to obtain periodic cell.
  std::unique_ptr<cl::CellMap> cell_map_;   //!< Related subcell indices to ranges of points.
};


}  // namespace cgrid


#endif  // CGRID_CELLGRID_H_

// vim: textwidth=90 et ts=2 sw=2
