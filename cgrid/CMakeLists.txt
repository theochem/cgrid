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

# Configure a header file to pass some of the CMake settings
# to the source code.
configure_file (
  ${CMAKE_CURRENT_SOURCE_DIR}/config.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/config.h
)

# Find cellcutoff include
find_path(CELL_INCLUDE_PATH cellcutoff/cell.h)

# Explain that the current directory has to be in the source path/
# This allows use to use include "cgrid/foo.h", which is consistent
# with the way external programs would compile against cgrid.
include_directories(${PROJECT_SOURCE_DIR} ${CELL_INCLUDE_PATH})

# Define source files
set(SOURCE_FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/cellgrid.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/scalarfns.cpp
)

# Define header files
set(HEADER_FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/cellgrid.h
  ${CMAKE_CURRENT_BINARY_DIR}/config.h
  ${CMAKE_CURRENT_SOURCE_DIR}/scalarfns.h
)

# Find the cellcutoff library. Let's be nice with manually compiled stuff in
# LD_LIBRARY_PATH...
find_library(CELLCUTOFF_LIBRARY cellcutoff
             HINTS ENV LD_LIBRARY_PATH)

# Define the shared library
add_library(cgrid SHARED ${SOURCE_FILES})
set_property(TARGET cgrid PROPERTY VERSION ${CGRID_VERSION})
set_property(TARGET cgrid PROPERTY SOVERSION ${CGRID_SOVERSION})
set_property(TARGET cgrid PROPERTY CXX_STANDARD 11)

# Link with cellcutoff
target_link_libraries(cgrid ${CELLCUTOFF_LIBRARY})

# Install include files
install(FILES ${HEADER_FILES} DESTINATION include/cgrid)

# Install library
install(TARGETS cgrid DESTINATION lib)

# Unit Tests with Google testing framework
add_subdirectory(tests)
