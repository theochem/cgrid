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

cmake_minimum_required(VERSION 3.0)
project(QCGrids CXX)

string(TOUPPER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE)
IF (CMAKE_BUILD_TYPE STREQUAL "DEBUG")
    #SET (CMAKE_VERBOSE_MAKEFILE "ON")
    IF (CMAKE_COMPILER_IS_GNUCXX)
        # When coding (debug mode), we want to know when there might be trouble ...
        SET (CMAKE_CXX_FLAGS "-g -O0 -Wall -Wconversion -Wextra -Werror -pedantic -coverage")
        SET (CMAKE_C_FLAGS "-g -O0 -Wall -Wconversion -Wextra -Werror -pedantic -coverage")
        SET (CMAKE_EXE_LINKER_FLAGS "-coverage")
        SET (QCGRIDS_DO_COVERAGE "ON")
    ENDIF ()
ENDIF ()


# Version information
set (QCGRIDS_MAJOR 0)
set (QCGRIDS_MINOR 0)
set (QCGRIDS_PATCH 0)

set (QCGRIDS_VERSION ${QCGRIDS_MAJOR}.${QCGRIDS_MINOR}.${QCGRIDS_PATCH})
set (QCGRIDS_SOVERSION ${QCGRIDS_MAJOR}.${QCGRIDS_MINOR})
set (PROJECT_BRIEF "A 3D domain decomposition library")

# Compilation of qcgrids according to qcgrids/CMakeLists.txt
add_subdirectory(qcgrids)
