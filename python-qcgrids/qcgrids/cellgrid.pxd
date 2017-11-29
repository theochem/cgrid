# -*- coding: utf-8 -*-
# QCGrids is a numerical integration library for quantum chemistry.
# Copyright (C) 2011-2017 The QCGrids developers
#
# This file is part of QCGrids.
#
# QCGrids is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# QCGrids is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>
# --

cimport cellcutoff.cell
from libcpp.vector cimport vector


cdef extern from "qcgrids/cellgrid.h" namespace "qcgrids":
    cdef cppclass CellgridPoint:
        CellgridPoint(const double* cart, const double weight, const int index);

    ctypedef double (*GridFunc) (const double*, const double, const void*);

    cdef cppclass Cellgrid:
        Cellgrid(const cellcutoff.cell.Cell& cell, double spacing);

        const vector[CellgridPoint]& grid_array();
        const cellcutoff.cell.Cell* cell();
        const cellcutoff.cell.Cell* subcell();
        const int* shape();

        void emplace_back_many(const double* cart, const double* weight, const size_t npoint);
        void sort();

        void iadd_cutoff(const double* center, const double cutoff,
            GridFunc grid_func, const void* extra_arg, double* output);
        double integrate_cutoff(const double* center, const double cutoff,
            GridFunc grid_func, const void* extra_arg, const double* factor);
