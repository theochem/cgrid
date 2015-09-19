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
#include <stdexcept>

#include <gtest/gtest.h>

#include <qcgrids/subgrid.h>

#include "common.h"


namespace qcg = qcgrids;


class SubgridTest : public ::testing::Test {
 public:
  virtual void SetUp() {
    const double center[3]{1.0, 0.0, -1.0};
    subgrid.reset(new qcg::Subgrid(center));
    double cart0[3]{4.0, 4.0, 0.0};
    double cart1[3]{-2.0, 4.0, 0.0};
    double cart2[3]{0.0, 1.0, 0.0};
    subgrid->emplace_back(cart0, 5.0, 0.2, 2);
    subgrid->emplace_back(cart1, 5.0, 0.3, 2);
    subgrid->emplace_back(cart2, sqrt(3.0), -0.5, 4);
  };

  std::unique_ptr<qcg::Subgrid> subgrid;
};


TEST_F(SubgridTest, iadd_super_example) {
  double fsuper[5]{1.0, 0.3, -0.2, 0.0, 0.8};
  double fsub[3]{2.0, -0.5, 0.2};
  subgrid->iadd_super(fsuper, fsub);
  EXPECT_NEAR(1.0, fsuper[0], EPS);
  EXPECT_NEAR(0.3, fsuper[1], EPS);
  EXPECT_NEAR(1.3, fsuper[2], EPS);
  EXPECT_NEAR(0.0, fsuper[3], EPS);
  EXPECT_NEAR(1.0, fsuper[4], EPS);
  EXPECT_NEAR(2.0, fsub[0], EPS);
  EXPECT_NEAR(-0.5, fsub[1], EPS);
  EXPECT_NEAR(0.2, fsub[2], EPS);
}


TEST_F(SubgridTest, take_sub_example) {
  double fsuper[5]{1.0, 0.3, -0.2, 0.0, 0.8};
  double fsub[3]{0.0, 0.0, 0.0};
  subgrid->take_sub(fsuper, fsub);
  EXPECT_NEAR(1.0, fsuper[0], EPS);
  EXPECT_NEAR(0.3, fsuper[1], EPS);
  EXPECT_NEAR(-0.2, fsuper[2], EPS);
  EXPECT_NEAR(0.0, fsuper[3], EPS);
  EXPECT_NEAR(0.8, fsuper[4], EPS);
  EXPECT_NEAR(-0.2, fsub[0], EPS);
  EXPECT_NEAR(-0.2, fsub[1], EPS);
  EXPECT_NEAR(0.8, fsub[2], EPS);
}


TEST_F(SubgridTest, integrate_example) {
  const double fsub0[3]{3.0, 2.0, -1.5};
  const double fsub1[6]{3.0, 4.0, -2.0, -5.0, -0.15, 0.1};
  const double *const factors[3]{fsub0, fsub1, fsub1+1};
  size_t strides[3]{1, 2, 2};
  double expected = 0.2*3.0*3.0*4.0 + 0.3*2.0*-2.0*-5.0 + -0.5*-1.5*-0.15*0.1;
  EXPECT_NEAR(expected, subgrid->integrate(factors, strides, 3), EPS);
}



// vim: textwidth=90 et ts=2 sw=2
