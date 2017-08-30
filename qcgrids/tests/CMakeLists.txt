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

# Make use of the GTest CMake module
find_package(GTest)

if(GTEST_FOUND)
    # Make sure the header files are found
    include_directories(${GTEST_INCLUDE_DIRS})

    # Find the cellcutoff library. Let's be nice with manually compiled stuff in
    # LD_LIBRARY_PATH...
    find_library(CELLCUTOFF_LIBRARY cellcutoff
                 HINTS ENV LD_LIBRARY_PATH)

    # Define test source files
    set(TEST_SOURCE_FILES
      ${CMAKE_CURRENT_SOURCE_DIR}/test_supergrid.cpp
      ${CMAKE_CURRENT_SOURCE_DIR}/common.cpp
    )

    # The testing executable
    add_executable(test_qcgrids ${TEST_SOURCE_FILES})

    # How the testing executable has to be linked...
    target_link_libraries(test_qcgrids qcgrids ${CELLCUTOFF_LIBRARY} ${GTEST_BOTH_LIBRARIES} gcov pthread)

    # Install the test in libexec
    install (TARGETS test_qcgrids DESTINATION libexec/qcgrids)

    # Use our own test target because ctest hides relevant output.
    if(QCGRIDS_DO_COVERAGE)
        add_custom_target(clean_coverage
                          COMMAND ${PROJECT_SOURCE_DIR}/tools/coverage_clean.sh
                          WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
                          COMMENT "Removing old coverage statistics")

        add_custom_target(check
                          COMMAND test_qcgrids
                          DEPENDS test_qcgrids clean_coverage
                          COMMENT "Running tests")

        add_custom_target(coverage
                          COMMAND ${PROJECT_SOURCE_DIR}/tools/coverage_analyze.py
                          DEPENDS check
                          WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
                          COMMENT "Coverage analysis")

    else()
        add_custom_target(check
                          COMMAND test_qcgrids
                          DEPENDS test_qcgrids
                          WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                          COMMENT "Running tests")


    endif()
endif()
