#!/usr/bin/env bash
./tools/cpplint.py --linelength=90 qcgrids/*.cpp qcgrids/*.h
./tools/cpplint.py --linelength=120 qcgrids/tests/*.cpp qcgrids/tests/*.h
