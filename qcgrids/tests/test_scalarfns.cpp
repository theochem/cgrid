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



#include <stdexcept>

#include <cmath>

#include "gtest/gtest.h"

#include "./common.h"
#include "qcgrids/scalarfns.h"



namespace qcg = qcgrids;


TEST(ExponentialTest, basics) {
  double exponent = 1.5;
  double x = 1.2;
  double y = 0.75971;
  qcg::Exponential a(exponent);
  EXPECT_NEAR(a.calc(x), exp(x*exponent), EPS);
  EXPECT_NEAR(a.calc_inv(a.calc(x)), x, EPS);
  EXPECT_NEAR(a.calc(a.calc_inv(y)), y, EPS);
  double output[3];
  a.calc(x, 2, output);
  EXPECT_NEAR(output[0], a.calc(x), EPS);
  EXPECT_NEAR(output[1], a.calc(x)*exponent, EPS);
  EXPECT_NEAR(output[2], a.calc(x)*exponent*exponent, EPS);
  a.calc_inv(y, 2, output);
  EXPECT_NEAR(output[0], a.calc_inv(y), EPS);
  EXPECT_NEAR(output[1], 1.0/(y*exponent), EPS);
  EXPECT_NEAR(output[2], -1.0/(y*y*exponent), EPS);
}


// vim: textwidth=90 et ts=2 sw=2
