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
    EXPECT_THROW(supergrid->iadd_cutoff(nullptr, 0.0, nullptr, nullptr, nullptr), std::logic_error);
    EXPECT_THROW(supergrid->integrate_cutoff(nullptr, 0.0, nullptr, nullptr, nullptr), std::logic_error);

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
    EXPECT_THROW(supergrid->iadd_cutoff(nullptr, 0.0, nullptr, nullptr, nullptr), std::logic_error);
    EXPECT_THROW(supergrid->integrate_cutoff(nullptr, 0.0, nullptr, nullptr, nullptr), std::logic_error);

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


double test_fn(const double* delta, const double d, const void* extra_arg) {
  const double& amplitude = reinterpret_cast<const double*>(extra_arg)[0];
  const double& exponent = reinterpret_cast<const double*>(extra_arg)[1];
  EXPECT_NEAR(d, vec3::norm(delta), EPS);
  return delta[0]*amplitude*exp(-exponent*d);
}


TEST_F(SupergridTest, iadd_integrate_cutoff_0) {
  size_t ncell_total = 0;
  double work[NPOINT];
  for (int irep = 0; irep < NREP; ++irep) {
    std::unique_ptr<qcg::Supergrid> supergrid(create_case_0(irep));
    ncell_total += supergrid->cell_map()->size();

    // Select a cutoff sphere
    double center[3];
    double cutoff;
    unsigned int seed = fill_random_double(irep + NREP, center, 3, -2, 2);
    seed = fill_random_double(seed, &cutoff, 1, 0.5, 2.5);

    // Set extra arguments
    double extra_arg[2];
    seed = fill_random_double(seed, extra_arg, 2, 1.0, 3.0);

    // Call iadd
    std::fill(work, work+NPOINT, 0.0);
    supergrid->iadd_cutoff(center, cutoff, test_fn, extra_arg, work);

    // Check every point in work, and computing integral
    double expected_integral = 0.0;
    size_t ipoint = 0;
    for (const qcg::SupergridPoint& point : supergrid->grid_array()) {
      double delta[3];
      vec3::delta(center, point.cart_, delta);
      double distance = vec3::norm(delta);
      if (distance < cutoff) {
        double expected = test_fn(delta, distance, extra_arg);
        EXPECT_NEAR(expected, work[ipoint], EPS);
        expected_integral += work[ipoint]*point.weight_;
      } else {
        EXPECT_EQ(0.0, work[ipoint]);
        // Set this work element to one, just to make sure it is not used in
        // integrate_cutoff (test below).
        work[ipoint] = 1.0;
      }
      ++ipoint;
    }

    // Compute integral in different ways
    double integral1 = supergrid->integrate_cutoff(center, cutoff, nullptr, nullptr, work);
    EXPECT_NEAR(expected_integral, integral1, EPS);
    double integral2 = supergrid->integrate_cutoff(center, cutoff, test_fn, extra_arg, nullptr);
    EXPECT_NEAR(expected_integral, integral2, EPS);
    std::fill(work, work+NPOINT, 2.0);
    double integral3 = supergrid->integrate_cutoff(center, cutoff, test_fn, extra_arg, work);
    EXPECT_NEAR(2*expected_integral, integral3, EPS);
  }

  // Sufficiency checks
  EXPECT_LT(3*NREP, ncell_total);
}


TEST_F(SupergridTest, iadd_integrate_cutoff_3) {
  size_t ncell_total = 0;
  size_t nnotvisited = 0;
  double vecs[9]{1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
  cl::Cell cell(vecs, 3);
  double work[NPOINT];
  for (int irep = 0; irep < NREP; ++irep) {
    std::unique_ptr<qcg::Supergrid> supergrid(create_case_3(irep, cell));
    ncell_total += supergrid->cell_map()->size();

    // Select a cutoff sphere
    double center[3];
    double cutoff;
    unsigned int seed = fill_random_double(irep + NREP, center, 3, -2, 2);
    seed = fill_random_double(seed, &cutoff, 1, 0.5, 2.5);

    // Set extra arguments
    double extra_arg[2];
    seed = fill_random_double(seed, extra_arg, 2, 1.0, 3.0);

    // Call iadd
    std::fill(work, work+NPOINT, 0.0);
    supergrid->iadd_cutoff(center, cutoff, test_fn, extra_arg, work);

    // Check every point in work, and computing integral
    double expected_integral = 0.0;
    size_t ipoint = 0;
    for (const qcg::SupergridPoint& point : supergrid->grid_array()) {
      double delta0[3];
      vec3::delta(center, point.cart_, delta0);
      int ranges_begin[3];
      int ranges_end[3];
      supergrid->cell()->ranges_cutoff(center, cutoff, ranges_begin, ranges_end);
      int icell[3];
      double expected = 0.0;
      bool visited = false;
      for (icell[0] = ranges_begin[0]; icell[0] < ranges_end[0]; ++icell[0]) {
        for (icell[1] = ranges_begin[1]; icell[1] < ranges_end[1]; ++icell[1]) {
          for (icell[2] = ranges_begin[2]; icell[2] < ranges_end[2]; ++icell[2]) {
            double delta[3];
            vec3::copy(delta0, delta);
            supergrid->cell()->iadd_vec(delta, icell);
            double distance = vec3::norm(delta);
            if (distance < cutoff) {
              expected += test_fn(delta, distance, extra_arg);
              visited = true;
            }
          }
        }
      }
      EXPECT_NEAR(expected, work[ipoint], EPS);
      expected_integral += work[ipoint]*point.weight_;
      if (!visited) {
        // Set this work element to one, just to make sure it is not used in
        // integrate_cutoff (test below).
        work[ipoint] = 1.0;
        ++nnotvisited;
      }
      ++ipoint;
    }

    // Compute integral in different ways
    /* The first one is not supposed to wrok as in the aperiodic case. The integrand is
       already wrapped periodically in work, and the following would again periodically
       repeast that wrapped result. This effectively "double wraps" the integrand, which
       does not make any sense.
    double integral1 = supergrid->integrate_cutoff(center, cutoff, nullptr, nullptr, work);
    EXPECT_NEAR(expected_integral, integral1, EPS);
    */
    double integral2 = supergrid->integrate_cutoff(center, cutoff, test_fn, extra_arg, nullptr);
    EXPECT_NEAR(expected_integral, integral2, EPS);
    std::fill(work, work+NPOINT, 2.0);
    double integral3 = supergrid->integrate_cutoff(center, cutoff, test_fn, extra_arg, work);
    EXPECT_NEAR(2*expected_integral, integral3, EPS*10);  // Order-of-operations differs
  }

  // Sufficiency checks
  EXPECT_LT(3*NREP, ncell_total);
  EXPECT_LT(0, nnotvisited);
}


TEST_F(SupergridTest, iadd_integrate_cutoff_errors) {
    std::unique_ptr<qcg::Supergrid> supergrid(new qcg::Supergrid(cl::Cell()));
    EXPECT_THROW(supergrid->iadd_cutoff(nullptr, 1.0, nullptr, nullptr, nullptr), std::logic_error);
    EXPECT_THROW(supergrid->integrate_cutoff(nullptr, 1.0, nullptr, nullptr, nullptr), std::logic_error);
}


// vim: textwidth=90 et ts=2 sw=2
