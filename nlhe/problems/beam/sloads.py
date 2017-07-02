#!/usr/bin/env python
# -*- coding: utf-8 -*-

import fenics as fe

from nlhe.problems.elproblem import SurfaceLoads

__author__ = "Alessio Nava"
__copyright__ = "2017 Alessio Nava"
__credits__ = ["Alessio Nava, the FEniCS Team"]
__license__ = "GPL"
__version__ = "3.0"
__maintainer__ = "Alessio Nava"
__email__ = "alessio.nava2@mail.polimi.it"
__status__ = "Development"

surface_loads = SurfaceLoads()

p = 500000.
q = p * 100.
surface_loads.add_surface_load(3, fe.Constant((0., p, 0.)))
surface_loads.add_surface_load(4, fe.Constant((0., 0., -q)))
