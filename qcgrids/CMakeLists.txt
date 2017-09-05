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

# Configure a header file to pass some of the CMake settings
# to the source code.
configure_file (
  ${CMAKE_CURRENT_SOURCE_DIR}/config.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/config.h
)

# Explain that the current directory has to be in the source path/
# This allows use to use include "qcgrids/foo.h", which is consistent
# with the way external programs would compile against qcgrids.
include_directories(${PROJECT_SOURCE_DIR})

# Define source files
set(SOURCE_FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/cellgrid.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/scalarfns.cpp
)

# Define header files
set(HEADER_FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/cellgrid.h
  ${CMAKE_CURRENT_BINARY_DIR}/config.h
  ${CMAKE_CURRENT_BINARY_DIR}/scalarfns.h
)

# Find the cellcutoff library. Let's be nice with manually compiled stuff in
# LD_LIBRARY_PATH...
find_library(CELLCUTOFF_LIBRARY cellcutoff
             HINTS ENV LD_LIBRARY_PATH)

# Define the shared library
add_library(qcgrids SHARED ${SOURCE_FILES})
set_property(TARGET qcgrids PROPERTY VERSION ${QCGRIDS_VERSION})
set_property(TARGET qcgrids PROPERTY SOVERSION ${QCGRIDS_SOVERSION})
set_property(TARGET qcgrids PROPERTY CXX_STANDARD 11)

# Link with cellcutoff
target_link_libraries(qcgrids ${CELLCUTOFF_LIBRARY})

# Install include files
install(FILES ${HEADER_FILES} DESTINATION include/qcgrids)

# Install library
install(TARGETS qcgrids DESTINATION lib)

# Unit Tests with Google testing framework
add_subdirectory(tests)
