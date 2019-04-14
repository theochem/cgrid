# -*- coding: utf-8 -*-
# CGrid is a library for for molecular numerical integration.
# Copyright (C) 2011-2019 The CGrid Development Team
#
# This file is part of CGrid.
#
# CGrid is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# CGrid is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>
# --


find_package(GTest REQUIRED)

# Define test source files
set(TEST_SOURCE_FILES
${CMAKE_CURRENT_SOURCE_DIR}/test_scalarfns.cpp
${CMAKE_CURRENT_SOURCE_DIR}/test_supergrid.cpp
${CMAKE_CURRENT_SOURCE_DIR}/common.cpp
)

# The testing executable
add_executable(test_cgrid ${TEST_SOURCE_FILES})
set_property(TARGET test_cgrid PROPERTY CXX_STANDARD 11)
set_property(TARGET test_cgrid PROPERTY include_directories ${GTEST_INCLUDE_DIRS})
# Google test uses many tricks that generate warnings, which cannot be fixed easily.
set_property(TARGET test_cgrid PROPERTY COMPILE_FLAGS -Wno-error)

target_link_libraries(test_cgrid cgrid ${CELLCUTOFF_LIBRARY} GTest::GTest GTest::Main pthread)

# Install the test in libexec
install (TARGETS test_cgrid DESTINATION libexec/cgrid)

# Use our own test target because ctest hides relevant output.
add_custom_target(test
                  COMMAND test_cgrid
                  DEPENDS test_cgrid
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                  COMMENT "Running tests")
