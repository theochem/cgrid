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

cmake_minimum_required(VERSION 3.1)
project(CGrid CXX)

# Version information
include(CMakeListsVersion.txt.in)
set(CGRID_VERSION ${GIT_TAG_VERSION})
set(CGRID_SOVERSION ${GIT_TAG_SOVERSION})
set(PROJECT_BRIEF "A library for for molecular numerical integration.")

# General compiler flags
set(CMAKE_CXX_FLAGS "-std=c++11 -Wall -Wconversion -Wextra")

# Additional compiler settings when doing a debug build
string(TOUPPER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE)
if(CMAKE_BUILD_TYPE STREQUAL "DEBUG")
    # These settings are only useful in dev mode.
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror -pedantic")
    if(CMAKE_COMPILER_IS_GNUCXX)
        # When coding (debug mode), we want to know when there might be trouble ...
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 -coverage")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_CXX_FLAGS} -coverage")
        set(CGRID_DO_COVERAGE "ON")
    endif()
endif()

# Compilation of cgrid according to cgrid/CMakeLists.txt
add_subdirectory(cgrid)

find_package(Doxygen)
if(DOXYGEN_FOUND)
    configure_file(${PROJECT_SOURCE_DIR}/doxygen.conf.in ${PROJECT_BINARY_DIR}/doxygen.conf @ONLY)
    add_custom_target(doc ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/doxygen.conf
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
    COMMENT "Generating API documentation with Doxygen" VERBATIM
)
endif(DOXYGEN_FOUND)

set(CPACK_SOURCE_GENERATOR "TGZ")
set(CPACK_SOURCE_PACKAGE_FILE_NAME "cgrid-${CGRID_VERSION}")
set(CPACK_SOURCE_IGNORE_FILES
  "/build/;/release/;/debug/;/.git/;/python-cgrid/;/HEADER/;${CPACK_SOURCE_IGNORE_FILES}")
include(CPack)
add_custom_target(sdist COMMAND ${CMAKE_MAKE_PROGRAM} package_source)
