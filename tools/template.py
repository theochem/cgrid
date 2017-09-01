#!/usr/bin/env python

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
