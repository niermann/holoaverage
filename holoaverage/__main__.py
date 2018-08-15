#!/usr/bin/env python3

# This file is part of holoaverage.
# Copyright (c) 2018 Tore Niermann
#
# holoaverage is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# holoaverage is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with holoaverage.  If not, see <http://www.gnu.org/licenses/>.

import sys
import os.path

# Insert this package directory into PYTHONPATH
# Clumsy but we assume if the user runs explicitly this package, he want to import from here and not from
# some other path (see also https://stackoverflow.com/questions/16981921/relative-imports-in-python-3)
PACKAGE_PARENT = '..'
SCRIPT_DIR = os.path.dirname(os.path.realpath(os.path.join(os.getcwd(), os.path.expanduser(__file__))))
sys.path.insert(0, os.path.normpath(os.path.join(SCRIPT_DIR, PACKAGE_PARENT)))

from holoaverage.main import main

if __name__ == '__main__':
    sys.exit(main(sys.argv))
