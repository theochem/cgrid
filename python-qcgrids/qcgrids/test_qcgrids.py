# -*- coding: utf-8 -*-
# QCGrids is a numerical integration library for quantum chemistry.
# Copyright (C) 2011-2015 Toon Verstraelen
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
#--


from qcgrids import *
import numpy as np

def test_subgrid_center():
    sg = Subgrid(np.array([1.0, -2.5, 3.0]))
    center = sg.center
    assert (center == np.array([1.0, -2.5, 3.0])).all()
    del sg
    assert (center == np.array([1.0, -2.5, 3.0])).all()


def test_subgrid_grid_array():
    sg = Subgrid(np.array([1.0, -2.5, 3.0]))
    sg.append(np.array([1.0, 2.0, 4.0]), np.sqrt(21.0), 0.3, 5)
    sg.append(np.array([3.0, 7.0, 9.0]), 132.0, 0.2, 6)
    sg.append(np.array([5.0, 8.0, 2.0]), 471.5714, 0.1, 7)
    assert len(sg.grid_array) == 3
    assert (sg.grid_array['cart'][0] == [1.0, 2.0, 4.0]).all()
    assert sg.grid_array['distance'][0] == np.sqrt(21.0)
    assert sg.grid_array['weight'][0] == 0.3
    assert sg.grid_array['index'][0] == 5
    print sg.grid_array['cart'][-1]
    assert (sg.grid_array['cart'][-1] == [5.0, 8.0, 2.0]).all()
    assert sg.grid_array['distance'][-1] == 471.5714
    assert sg.grid_array['weight'][-1] == 0.1
    assert sg.grid_array['index'][-1] == 7


def test_subgrid_grid_parts():
    sg = Subgrid(np.array([1.0, -2.5, 3.0]))
    sg.append(np.array([1.0, 2.0, 4.0]), np.sqrt(21.0), 0.3, 5)
    sg.append(np.array([3.0, 7.0, 9.0]), 132.0, 0.2, 6)
    sg.append(np.array([5.0, 8.0, 2.0]), 471.5714, 0.1, 7)
    assert len(sg.points) == 3
    assert len(sg.distances) == 3
    assert len(sg.weights) == 3
    assert len(sg.indices) == 3
    assert (sg.points[0] == [1.0, 2.0, 4.0]).all()
    assert sg.distances[0] == np.sqrt(21.0)
    assert sg.weights[0] == 0.3
    assert sg.indices[0] == 5
    assert (sg.points[-1] == [5.0, 8.0, 2.0]).all()
    assert sg.distances[-1] == 471.5714
    assert sg.weights[-1] == 0.1
    assert sg.indices[-1] == 7
