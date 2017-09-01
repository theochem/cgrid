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


#include <limits>
#include <random>
#include <stdexcept>

#include <gtest/gtest.h>

#include "common.h"


//! Internal helper that just makes a useful random seed for a random generates.
unsigned int get_next_seed(std::minstd_rand gen) {
  std::uniform_int_distribution<unsigned int>
      dis_seed(0, std::numeric_limits<unsigned int>::max());
  return dis_seed(gen);
}


unsigned int fill_random_double(const unsigned int seed, double* array, const int size,
    const double low, const double high) {
  // Parameter check
  if (size <= 0)
    throw std::domain_error("Array size must be strictly positive.");

  // Fill the array with random data for given seed.
  std::minstd_rand gen(seed);
  std::uniform_real_distribution<double> dis(low, high);
  for (int i=0; i < size; ++i)
    array[i] = dis(gen);

  // Generate a different seed for the next call
  return get_next_seed(gen);
}

/*
unsigned int fill_random_int(const unsigned int seed, int* array, const int size,
    const int begin, const int end) {
  // Parameter check
  if (size <= 0)
    throw std::domain_error("Array size must be strictly positive.");
  if (begin > end)
    throw std::domain_error("Begin cannot be larger than end.");

  // Fill the array with random data for given seed.
  std::minstd_rand gen(seed);
  std::uniform_int_distribution<int> dis(begin, end);
  for (int i=0; i < size; ++i)
    array[i] = dis(gen);

  // Generate a different seed for the next call
  return get_next_seed(gen);
}
*/

TEST(CommonTest, domain) {
  EXPECT_THROW(fill_random_double(0, nullptr, 0, 0.0, 1.0), std::domain_error);
  EXPECT_THROW(fill_random_double(0, nullptr, -1, 0.0, 1.0), std::domain_error);
  // EXPECT_THROW(fill_random_int(0, nullptr, 0, 0, 1), std::domain_error);
  // EXPECT_THROW(fill_random_int(0, nullptr, -1, 0, 1), std::domain_error);
  // EXPECT_THROW(fill_random_int(0, nullptr, 1, 1, 0), std::domain_error);
}

// vim: textwidth=90 et ts=2 sw=2
