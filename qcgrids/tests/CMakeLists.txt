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


# Define test source files
set(TEST_SOURCE_FILES
${CMAKE_CURRENT_SOURCE_DIR}/test_scalarfns.cpp
${CMAKE_CURRENT_SOURCE_DIR}/test_supergrid.cpp
${CMAKE_CURRENT_SOURCE_DIR}/common.cpp
)

# The testing executable
add_executable(test_qcgrids ${TEST_SOURCE_FILES})
set_property(TARGET test_qcgrids PROPERTY CXX_STANDARD 11)
# Google test uses many tricks that generate warnings, which cannot be fixed easily.
set_property(TARGET test_qcgrids PROPERTY COMPILE_FLAGS -Wno-error)

# How the testing executable has to be linked...
if(UNIX AND NOT APPLE)
  set(GCOV gcov)
endif()

target_link_libraries(test_qcgrids qcgrids ${CELLCUTOFF_LIBRARY} gtest gtest_main ${GCOV} pthread)

# Install the test in libexec
install (TARGETS test_qcgrids DESTINATION libexec/qcgrids)

# Use our own test target because ctest hides relevant output.
add_custom_target(check
                  COMMAND test_qcgrids
                  DEPENDS test_qcgrids
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                  COMMENT "Running tests")
