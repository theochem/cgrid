#!/usr/bin/env bash

if [[ ${BUILD_TYPE} ]]; then
  BTYPE=${BUILD_TYPE}
else
  BTYPE="release"
fi

mkdir build
cd build
cmake \
  -DCMAKE_INSTALL_PREFIX=${PREFIX} \
  -DCMAKE_BUILD_TYPE=${BTYPE} \
  ..
make install
