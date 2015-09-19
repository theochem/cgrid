#!/usr/bin/env bash

rm -fvr gcov
find -type f | grep "gcda$" | xargs rm -vf
