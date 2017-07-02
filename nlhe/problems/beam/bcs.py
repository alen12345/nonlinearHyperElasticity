#!/usr/bin/env python
# -*- coding: utf-8 -*-

import fenics as fe

from nlhe.problems.elproblem import BoundaryConditions

__author__ = "Alessio Nava"
__copyright__ = "2017 Alessio Nava"
__credits__ = ["Alessio Nava, the FEniCS Team"]
__license__ = "GPL"
__version__ = "3.0"
__maintainer__ = "Alessio Nava"
__email__ = "alessio.nava2@mail.polimi.it"
__status__ = "Development"

boundary_conditions = BoundaryConditions()

boundary_conditions.add_boundary_condition(2, fe.Constant((0., 0., 0.)))
