#!/bin/bash
echo "Cleaning code in ${PWD} and subdirectories."
# Remove files that "contaminate" the source tree. The cmake build dir should be included
# here.
rm -vf python-qcgrids/qcgrids/*.pyc
rm -vf python-qcgrids/qcgrids/qcgrids.cpp
rm -vf python-qcgrids/MANIFEST
rm -vfr python-qcgrids/build
rm -vfr python-qcgrids/dist
# split output of find at newlines.
IFS=$'\n'
# send all relevant files to the code cleaner
find qcgrids *.* tools python-qcgrids | \
    egrep "(\.cpp$)|(\.h$)|(\.in$)|(\.sh$)|(\.py$)|(\.pyx$)|(\.pxd$)|(\.txt$)|(\.conf$)|(.gitignore$)" | \
    xargs ./tools/codecleaner.py
