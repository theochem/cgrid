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
#include <memory>
#include <iostream>

#include <cmath>

#include "gtest/gtest.h"

#include "./common.h"
#include "./ridders.h"
#include "qcgrids/scalarfns.h"



namespace qcg = qcgrids;


/** @brief
      A scalar function, with an x, y and f(x) value, which should be used for
      parameterized test cases.
*/

class SFTestParams {
 public:
  //! No default constructor allowed.
  SFTestParams() = delete;
  /** @brief
        Create a SFTestParams object
  */
  SFTestParams(qcg::ScalarFunction* sfn, double x, double y, double val) :
      sfn(sfn), x(x), y(y), val(val) {}

  const qcg::ScalarFunction* sfn;
  const double x;
  const double y;
  const double val;
};

class ScalarFunctionTest : public ::testing::TestWithParam<SFTestParams> {};


TEST_P(ScalarFunctionTest, basics) {
  const qcg::ScalarFunction* sfn(GetParam().sfn);
  const double x(GetParam().x);
  const double y(GetParam().y);
  EXPECT_NEAR(sfn->value(x), GetParam().val, EPS);
  {
    // Test consistency between general and specialized functions.
    EXPECT_NEAR(sfn->qcg::ScalarFunction::value(x), sfn->value(x), EPS);
    EXPECT_NEAR(sfn->qcg::ScalarFunction::deriv(x), sfn->deriv(x), EPS);
    EXPECT_NEAR(sfn->qcg::ScalarFunction::deriv2(x), sfn->deriv2(x), EPS);
    // Test consistency between array output and individual functions.
    double output[3];
    sfn->calc(x, 2, output);
    EXPECT_NEAR(output[0], sfn->value(x), EPS);
    EXPECT_NEAR(output[1], sfn->deriv(x), EPS);
    EXPECT_NEAR(output[2], sfn->deriv2(x), EPS);
    // Test derivative with finite-difference approximation
    double deriv = diff_ridders(sfn, &qcg::ScalarFunction::value, x, 0.1);
    EXPECT_NEAR(output[1], deriv, EPS);
    double deriv2 = diff_ridders(sfn, &qcg::ScalarFunction::deriv, x, 0.1);
    EXPECT_NEAR(output[2], deriv2, EPS);
    // Test some exceptions
    EXPECT_THROW(sfn->calc(x, -1, output), std::domain_error);
    EXPECT_THROW(sfn->calc(x, 3, output), std::domain_error);
  }
  if (sfn->invertible()) {
    // Test consistency of function and its inverse.
    EXPECT_NEAR(sfn->value_inv(sfn->value(x)), x, EPS);
    EXPECT_NEAR(sfn->value(sfn->value_inv(y)), y, EPS);
    // Test consistency between general and specialized functions.
    EXPECT_NEAR(sfn->qcg::ScalarFunction::value_inv(x), sfn->value_inv(x), EPS);
    EXPECT_NEAR(sfn->qcg::ScalarFunction::deriv_inv(x), sfn->deriv_inv(x), EPS);
    EXPECT_NEAR(sfn->qcg::ScalarFunction::deriv2_inv(x), sfn->deriv2_inv(x), EPS);
    // Test consistency between array output and individual functions.
    double output[3];
    sfn->calc_inv(y, 2, output);
    EXPECT_NEAR(output[0], sfn->value_inv(y), EPS);
    EXPECT_NEAR(output[1], sfn->deriv_inv(y), EPS);
    EXPECT_NEAR(output[2], sfn->deriv2_inv(y), EPS);
    // Test derivative with finite-difference approximation
    double deriv_inv = diff_ridders(sfn, &qcg::ScalarFunction::value_inv, y, 0.1);
    EXPECT_NEAR(output[1], deriv_inv, EPS);
    double deriv2_inv = diff_ridders(sfn, &qcg::ScalarFunction::deriv_inv, y, 0.1);
    EXPECT_NEAR(output[2], deriv2_inv, EPS);
    // Test some exceptions
    EXPECT_THROW(sfn->calc_inv(x, -1, output), std::domain_error);
    EXPECT_THROW(sfn->calc_inv(x, 3, output), std::domain_error);
  } else {
    EXPECT_THROW(sfn->value_inv(y), std::logic_error);
  }
}



qcg::Exp exp_case1(0.5, 1.5);
qcg::Exp exp_case2(0.1, 2.5);
qcg::Ln ln_case1(0.5, 0.2);
qcg::Ln ln_case2(0.1, 1.2);
qcg::Linear linear_case(0.4, 1.2);
qcg::Identity identity_case;
qcg::Constant constant_case(0.88);
qcg::Power power_case1(0.3, 2.3);
qcg::Power power_case2(1.1, -2);
qcg::Power power_case3(0.7, 1.0);
qcg::Rational rational_case(20.0, 33.0);
const double spline_input1[3] = {0.5, 1.0, 2.3};
const double spline_input2[3] = {-0.1, 1.2, 2.0};
qcg::UniformCubicSpline spline3(3, spline_input1, spline_input2);


INSTANTIATE_TEST_CASE_P(Examples, ScalarFunctionTest, ::testing::Values(
  SFTestParams(&exp_case1, 1.2, 0.75971, 0.5*exp(1.5*1.2)),
  SFTestParams(&exp_case2, 0.2, 0.8, 0.1*exp(2.5*0.2)),
  SFTestParams(&ln_case1, 1.2, 0.75971, 0.5*log(0.2*1.2)),
  SFTestParams(&ln_case2, 0.2, 0.6, 0.1*log(1.2*0.2)),
  SFTestParams(&linear_case, 0.2, 0.6, 0.4*0.2 + 1.2),
  SFTestParams(&identity_case, 0.2, 0.6, 0.2),
  SFTestParams(&constant_case, 0.2, 0.6, 0.88),
  SFTestParams(&power_case1, 0.5, 0.6, 0.3*pow(0.5, 2.3)),
  SFTestParams(&power_case2, 0.7, 1.3, 1.1*pow(0.7, -2)),
  SFTestParams(&power_case3, 0.5, 0.9, 0.7*0.5),
  SFTestParams(&rational_case, 11.0, 1.3, 20.0*11.0/(1 - 11.0/33.0)),
  SFTestParams(&spline3, 0.2, 0.3, 0.5008)
));


TEST(ScalarFunctionTest, corner_cases) {
  EXPECT_NEAR(power_case1.value(0.0), 0.0, EPS);
  EXPECT_NEAR(power_case1.value_inv(0.0), 0.0, EPS);
  EXPECT_NEAR(power_case1.deriv(0.0), 0.0, EPS);
  EXPECT_NEAR(power_case1.deriv2(0.0), 0.0, EPS);
}


TEST(ScalarFunctionTest, exceptions) {
  EXPECT_THROW(qcg::Exp(0.0, 1.0), std::domain_error);
  EXPECT_THROW(qcg::Exp(1.0, 0.0), std::domain_error);
  EXPECT_THROW(qcg::Ln(0.0, 1.0), std::domain_error);
  EXPECT_THROW(qcg::Ln(1.0, 0.0), std::domain_error);
  EXPECT_THROW(qcg::Linear(0.0, 1.0), std::domain_error);
  EXPECT_THROW(qcg::Power(0.0, 1.0), std::domain_error);
  EXPECT_THROW(qcg::Power(1.0, 0.0), std::domain_error);
  EXPECT_THROW(qcg::Rational(1.0, 0.0), std::domain_error);
  EXPECT_THROW(qcg::Rational(0.0, 1.0), std::domain_error);
  EXPECT_THROW(qcg::UniformCubicSpline(0), std::domain_error);
  EXPECT_THROW(qcg::UniformCubicSpline(1), std::domain_error);
}

TEST(ScalarFunctionTest, uniform_cubic_spline_basics) {
  EXPECT_EQ(spline3.npoint(), 3);
  EXPECT_EQ(spline3.left(), 0.0);
  EXPECT_EQ(spline3.right(), 2.0);
  EXPECT_EQ(spline3.x(0), 0.0);
  EXPECT_EQ(spline3.x(1), 1.0);
  EXPECT_EQ(spline3.x(2), 2.0);
  EXPECT_EQ(spline3.values()[0], 0.5);
  EXPECT_EQ(spline3.values()[1], 1.0);
  EXPECT_EQ(spline3.values()[2], 2.3);
  EXPECT_EQ(spline3.derivs()[0], -0.1);
  EXPECT_EQ(spline3.derivs()[1], 1.2);
  EXPECT_EQ(spline3.derivs()[2], 2.0);
}

// vim: textwidth=90 et ts=2 sw=2
