#!/usr/bin/env python
# -*- coding: utf-8 -*-
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
#!/usr/bin/env python

import os, subprocess, shutil, glob

# Create a clean working dir
if os.path.exists('gcov'):
    shutil.rmtree('gcov')
os.mkdir('gcov')
os.chdir('gcov')

# Run gcov on all relevant files
for root, dirs, fns in os.walk('../'):
    for fn in fns:
        if fn.endswith('.gcno'):
             fn_o = '%s/%s.o' % (root, fn[:-5])
             #print 'Coverage %s' % fn_o
             subprocess.call('gcov %s > /dev/null' % fn_o, shell=True)

# Sometimes gcov signals non-executed lines that are not relevant
def is_relevant(line):
    #print line[16:].strip()
    if line[16:].strip().startswith('class'):
        return False
    if line[16:].split() == ['do', '{']:
        return False
    return True

# Loop over all gcov files and extract relevant information
coverage_results = {}
for fn in os.listdir('.'):
    info_lines = []
    path = ''
    if fn.endswith('.gcov'):
        with open(fn) as f:
            first = f.next()[23:-1]
            if first.startswith('/usr') or first.endswith('.h'):
                continue
            path = first[23:-1]
            for line in f:
                if line.startswith('    #####:') and is_relevant(line):
                    info_lines.append(line[10:-1])
    dn, bn = os.path.split(path)
    l = coverage_results.setdefault(dn, [])
    l.append((bn, info_lines, fn))

# Nicely format things on screen
green='\033[32m'
endc = '\033[0m'
yellow = '\033[93m'
red = '\033[91m'
print 'Number of lines not covered in unit tests:'
for dn, bns in sorted(coverage_results.iteritems()):
    print '   ', dn
    for bn, info_lines, fn_gcov in sorted(bns):
        if len(info_lines) == 0:
            print '%s%8i  %s%s' % (green, len(info_lines), bn, endc)
        else:
            print '%s%8i  %s gcov/%s%s' % (red, len(info_lines), bn.ljust(30), fn_gcov, endc)
            for info_line in info_lines:
                print '         '+yellow+info_line+endc
