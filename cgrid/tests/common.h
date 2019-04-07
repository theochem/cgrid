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


#ifndef CGRID_TESTS_COMMON_H_
#define CGRID_TESTS_COMMON_H_


// Number of repetitions on a test with random data
#define NREP 100
// Number of random Cartesian points to sample in a test
#define NPOINT 1000
// Allowed tolerance when comparing doubles (EPS stands for epsilon.)
#define EPS 1e-10


/* Some notes:

    - Usually, the output argument is last. The exception in this module is the size
      of the output argument, which comes after the actual output argument.

*/


//! Fills an array of doubles with random numbers in range ]-0.5*scale, 0.5*scale]
unsigned int fill_random_double(const unsigned int seed, double* array, const int size,
    const double low = -0.5, const double high = 0.5);

/*
//! Fills an array of int with random numbers in range [-range, range]
unsigned int fill_random_int(const unsigned int seed, int* array, const int size,
    const int begin, const int end);
*/

#endif  // CGRID_TESTS_COMMON_H_

// vim: textwidth=90 et ts=2 sw=2
