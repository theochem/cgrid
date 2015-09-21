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


#include <cmath>
#include <memory>
#include <stdexcept>
#include <vector>

#include <celllists/cell.h>
#include <celllists/vec3.h>
#include <gtest/gtest.h>

#include <qcgrids/supergrid.h>
#include <qcgrids/subgrid.h>

#include "common.h"


namespace qcg = qcgrids;
namespace vec3 = celllists::vec3;


// Fixtures
// ========

class SupergridTest : public ::testing::Test {
 public:
  qcg::Supergrid* create_case_0(unsigned int seed) {
    // Construct the grid
    cl::Cell cell;
    qcg::Supergrid* supergrid = new qcg::Supergrid(cell);

    // Add points
    for (int ipoint = 0; ipoint < NPOINT; ipoint++) {
      double cart[3];
      double weight;
      seed = fill_random_double(seed, cart, 3, -5, 5);
      seed = fill_random_double(seed, &weight, 1, 1, 2);
      supergrid->emplace_back(cart, weight);
      EXPECT_EQ(cart[0], supergrid->grid_array()[ipoint].cart_[0]);
      EXPECT_EQ(cart[1], supergrid->grid_array()[ipoint].cart_[1]);
      EXPECT_EQ(cart[2], supergrid->grid_array()[ipoint].cart_[2]);
      EXPECT_EQ(weight, supergrid->grid_array()[ipoint].weight_);
      EXPECT_EQ(0, supergrid->grid_array()[ipoint].icell_[0]);
      EXPECT_EQ(0, supergrid->grid_array()[ipoint].icell_[1]);
      EXPECT_EQ(0, supergrid->grid_array()[ipoint].icell_[2]);
      EXPECT_EQ(ipoint, supergrid->grid_array()[ipoint].index_);
    }

    // Basic checks on data members
    EXPECT_EQ(NPOINT, supergrid->grid_array().size());
    EXPECT_TRUE(std::isnan(supergrid->cell()->volume()));
    EXPECT_NEAR(1.0, supergrid->subcell()->volume(), EPS);
    EXPECT_EQ(0, supergrid->shape()[0]);
    EXPECT_EQ(0, supergrid->shape()[1]);
    EXPECT_EQ(0, supergrid->shape()[2]);

    // Check for the right exception (sort isn't called yet)
    EXPECT_THROW(supergrid->cell_map(), std::logic_error);
    EXPECT_THROW(supergrid->create_subgrid(nullptr, 1.0), std::logic_error);

    // Sort the grid_array, assign icell and make cell_map
    supergrid->sort();

    // Test if points are sorted
    const std::vector<qcg::SupergridPoint>& grid_array(supergrid->grid_array());
    size_t ncell = 0;
    for (size_t ipoint = 1; ipoint < NPOINT; ipoint++) {
      EXPECT_LE(grid_array[ipoint - 1].icell_[0], grid_array[ipoint].icell_[0]);
      if (grid_array[ipoint - 1].icell_[0] == grid_array[ipoint].icell_[0]) {
        EXPECT_LE(grid_array[ipoint - 1].icell_[1], grid_array[ipoint].icell_[1]);
        if (grid_array[ipoint - 1].icell_[1] == grid_array[ipoint].icell_[1]) {
          EXPECT_LE(grid_array[ipoint - 1].icell_[2], grid_array[ipoint].icell_[2]);
        }
      }
      if ((grid_array[ipoint - 1].icell_[0] != grid_array[ipoint].icell_[0]) ||
          (grid_array[ipoint - 1].icell_[1] != grid_array[ipoint].icell_[1]) ||
          (grid_array[ipoint - 1].icell_[2] != grid_array[ipoint].icell_[2]))
        ++ncell;
    }
    ++ncell;

    // Basic test of the cell_map
    EXPECT_EQ(ncell, supergrid->cell_map()->size());

    return supergrid;
  }

  qcg::Supergrid* create_case_3(unsigned int seed, const cl::Cell& cell) {
    // Construct the grid
    qcg::Supergrid* supergrid = new qcg::Supergrid(cell, 0.55);

    // Add points
    double all_cart[3*NPOINT];
    double all_weight[NPOINT];
    seed = fill_random_double(seed, all_cart, 3*NPOINT, -5, 5);
    seed = fill_random_double(seed, all_weight, NPOINT, 1, 2);
    supergrid->emplace_back_many(all_cart, all_weight, NPOINT);
    for (int ipoint = 0; ipoint < NPOINT; ipoint++) {
      EXPECT_EQ(all_cart[3*ipoint], supergrid->grid_array()[ipoint].cart_[0]);
      EXPECT_EQ(all_cart[3*ipoint + 1], supergrid->grid_array()[ipoint].cart_[1]);
      EXPECT_EQ(all_cart[3*ipoint + 2], supergrid->grid_array()[ipoint].cart_[2]);
      EXPECT_EQ(all_weight[ipoint], supergrid->grid_array()[ipoint].weight_);
      EXPECT_EQ(0, supergrid->grid_array()[ipoint].icell_[0]);
      EXPECT_EQ(0, supergrid->grid_array()[ipoint].icell_[1]);
      EXPECT_EQ(0, supergrid->grid_array()[ipoint].icell_[2]);
      EXPECT_EQ(ipoint, supergrid->grid_array()[ipoint].index_);
    }

    // Basic checks on data members
    EXPECT_EQ(NPOINT, supergrid->grid_array().size());
    EXPECT_NEAR(1.0, supergrid->cell()->volume(), EPS);
    EXPECT_NEAR(0.125, supergrid->subcell()->volume(), EPS);
    EXPECT_EQ(2, supergrid->shape()[0]);
    EXPECT_EQ(2, supergrid->shape()[1]);
    EXPECT_EQ(2, supergrid->shape()[2]);

    // Check for the right exception (sort isn't called yet)
    EXPECT_THROW(supergrid->cell_map(), std::logic_error);
    EXPECT_THROW(supergrid->create_subgrid(nullptr, 1.0), std::logic_error);

    // Sort the grid_array, assign icell and make cell_map
    supergrid->sort();

    // Test if points are sorted
    const std::vector<qcg::SupergridPoint>& grid_array(supergrid->grid_array());
    size_t ncell = 0;
    for (size_t ipoint = 1; ipoint < NPOINT; ipoint++) {
      EXPECT_LE(grid_array[ipoint - 1].icell_[0], grid_array[ipoint].icell_[0]);
      if (grid_array[ipoint - 1].icell_[0] == grid_array[ipoint].icell_[0]) {
        EXPECT_LE(grid_array[ipoint - 1].icell_[1], grid_array[ipoint].icell_[1]);
        if (grid_array[ipoint - 1].icell_[1] == grid_array[ipoint].icell_[1]) {
          EXPECT_LE(grid_array[ipoint - 1].icell_[2], grid_array[ipoint].icell_[2]);
        }
      }
      if ((grid_array[ipoint - 1].icell_[0] != grid_array[ipoint].icell_[0]) ||
          (grid_array[ipoint - 1].icell_[1] != grid_array[ipoint].icell_[1]) ||
          (grid_array[ipoint - 1].icell_[2] != grid_array[ipoint].icell_[2]))
        ++ncell;
    }
    ++ncell;

    // Basic test of the cell_map
    EXPECT_EQ(ncell, supergrid->cell_map()->size());

    return supergrid;
  }
};


TEST_F(SupergridTest, subgrid_example1) {
  size_t ncell_total = 0;
  size_t npoint_inside = 0;
  for (int irep = 0; irep < NREP; ++irep) {
    std::unique_ptr<qcg::Supergrid> supergrid(create_case_0(irep));
    ncell_total += supergrid->cell_map()->size();

    // Make a subgrid
    double center[3];
    double cutoff;
    unsigned int seed = fill_random_double(irep + NREP, center, 3, -2, 2);
    fill_random_double(seed, &cutoff, 1, 3, 15);
    std::unique_ptr<qcg::Subgrid> subgrid(supergrid->create_subgrid(center, cutoff));

    // Check basics of subgrid
    EXPECT_EQ(center[0], subgrid->center()[0]);
    EXPECT_EQ(center[1], subgrid->center()[1]);
    EXPECT_EQ(center[2], subgrid->center()[2]);

    // Compare all the points in the subgrid with the supergrid
    for (const qcg::SubgridPoint& point : subgrid->grid_array()) {
      EXPECT_GT(supergrid->grid_array().size(), point.index_);
      // Cartesian coordinates must be the same because there are no periodic boundary
      // conditions.
      EXPECT_EQ(supergrid->grid_array()[point.index_].cart_[0], point.cart_[0]);
      EXPECT_EQ(supergrid->grid_array()[point.index_].cart_[1], point.cart_[1]);
      EXPECT_EQ(supergrid->grid_array()[point.index_].cart_[2], point.cart_[2]);
      // Weights should always be the same
      EXPECT_EQ(supergrid->grid_array()[point.index_].weight_, point.weight_);
      // Check correctness of distance
      EXPECT_GE(cutoff, point.distance_);
      EXPECT_NEAR(vec3::distance(center, point.cart_), point.distance_, EPS);
      ++npoint_inside;
    }

    // Go over all the points in the supergrid and check that each point within the cutoff
    // is included in the subgrid.
    int isuper = 0;
    for (const qcg::SupergridPoint point0 : supergrid->grid_array()) {
      const double d(vec3::distance(center, point0.cart_));
      if (d <= cutoff) {
        bool found(false);
        for (const qcg::SubgridPoint point1 : subgrid->grid_array()) {
          if (point1.index_ == isuper) {
            found = true;
            break;
          }
        }
        EXPECT_EQ(true, found);
      }
      ++isuper;
    }
  }

  // Sufficiency checks
  EXPECT_LT(3*NREP, ncell_total);
  EXPECT_LT((NREP*NPOINT)/10, npoint_inside);
}


TEST_F(SupergridTest, subgrid_example2) {
  size_t ncell_total = 0;
  size_t npoint_inside = 0;
  double vecs[9]{1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
  cl::Cell cell(vecs, 3);
  for (int irep = 0; irep < NREP; ++irep) {
    std::unique_ptr<qcg::Supergrid> supergrid(create_case_3(irep, cell));
    ncell_total += supergrid->cell_map()->size();

    // Make a subgrid
    double center[3];
    double cutoff;
    unsigned int seed = fill_random_double(irep + NREP, center, 3, -2, 2);
    fill_random_double(seed, &cutoff, 1, 0.5, 2.5);
    std::unique_ptr<qcg::Subgrid> subgrid(supergrid->create_subgrid(center, cutoff));

    // Check basics of subgrid
    EXPECT_EQ(center[0], subgrid->center()[0]);
    EXPECT_EQ(center[1], subgrid->center()[1]);
    EXPECT_EQ(center[2], subgrid->center()[2]);

    // Compare all the points in the subgrid with the supergrid
    for (const qcg::SubgridPoint& point : subgrid->grid_array()) {
      EXPECT_GT(supergrid->grid_array().size(), point.index_);
      // The cartesian coordinates may differ up to a linear combination of cell vectors.
      // Hence, only tesing wrapped relative fractional coordinates
      double delta[3];
      vec3::delta(supergrid->grid_array()[point.index_].cart_, point.cart_, delta);
      cell.iwrap_mic(delta);
      EXPECT_NEAR(0.0, delta[0], EPS);
      EXPECT_NEAR(0.0, delta[1], EPS);
      EXPECT_NEAR(0.0, delta[2], EPS);
      // Weights should always be the same
      EXPECT_EQ(supergrid->grid_array()[point.index_].weight_, point.weight_);
      // Check correctness of distance
      EXPECT_GE(cutoff, point.distance_);
      EXPECT_NEAR(vec3::distance(center, point.cart_), point.distance_, EPS);
      ++npoint_inside;
    }
  }

  // Sufficiency checks
  EXPECT_LT(3*NREP, ncell_total);
  EXPECT_LT((NREP*NPOINT)/10, npoint_inside);
}


// vim: textwidth=90 et ts=2 sw=2
