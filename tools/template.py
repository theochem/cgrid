#!/usr/bin/env python
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

from __future__ import print_function

from glob import glob
from string import Template

import yaml


def main():
    for template_fn in glob('*-template.*') + glob('.*-template.*'):
        if template_fn.endswith('.vars.yml'):
            continue
        print("Processing {}".format(template_fn))
        # Load the template
        with open(template_fn) as f:
            template = Template(f.read())
        # Load the variables
        with open(template_fn + '.vars.yml') as f:
            variables = yaml.load(f)
        with open('../' + template_fn.replace('-template', ''), "w") as f:
            f.write(template.safe_substitute(**variables))

if __name__ == '__main__':
    main()
