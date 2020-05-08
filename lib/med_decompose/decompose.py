#!/usr/bin/env python3
# -*- coding: utf-8

# Copyright (C) 2010 - 2012, A. Murat Eren
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# Please read the COPYING file.

import os
import sys

from Oligotyping.lib.decomposer import Decomposer
from Oligotyping.utils.utils import ConfigError
from Oligotyping.utils import parsers


# med_out_dir = '/Users/humebc/phylogenetic_software/med_hume'
# path_to_file = '/Users/humebc/phylogenetic_software/med_hume/test_fasta.fasta'

# ['-M', '4', '--skip-gexf-files', '--skip-gen-figures', '--skip-gen-html', '--skip-check-input', '-o', med_out_dir,
#  path_to_file]
parser = parsers.decomposer()
decomposer = Decomposer(parser.parse_args())

if __name__ == '__main__':

    try:
        decomposer.decompose()
    except ConfigError as e:
        print(e)
        sys.exit(-1)

