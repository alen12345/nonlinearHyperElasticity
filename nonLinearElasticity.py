#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
import fenics as fe
import solver.nlhe as nlhe

__author__ = "Alessio Nava"
__copyright__ = "2017 Alessio Nava"
__credits__ = ["Alessio Nava, the FEniCS Team"]
__license__ = "GPL"
__version__ = "3.0"
__maintainer__ = "Alessio Nava"
__email__ = "alessio.nava2@mail.polimi.it"
__status__ = "Development"

# Mesh name
mesh_name = "beam"

# Material properties
materials = {1: [210e9, 0.33],
             }

# Define the boundary conditions
boundary_conditions = {2: fe.Constant((0., 0., 0.)),
                       }

# Load definitions
p = 500000.
q = p * 100.
surface_loads = {3: fe.Constant((0., p, 0.)),
                 4: fe.Constant((0., 0., -q))
                 }

body_forces = {1: fe.Constant((0., 7800*9.81, 0.)),
               }


# Stored strain energy density for a generic model
def strainDensityFunction(E, Ey, nu):
    mu = Ey / (2.*(1+nu))
    lambda_ = Ey*nu / ((1+nu)*(1-2*nu))
    return lambda_/2.*(fe.tr(E))**2. + mu*fe.tr(E*E)

# Invoke the solver
nlhe.nonlinearHyperElasticitySolver(mesh_name, materials, boundary_conditions,
                                    surface_loads, body_forces,
                                    strainDensityFunction)