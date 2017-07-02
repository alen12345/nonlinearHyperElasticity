#!/usr/bin/env python
# -*- coding: utf-8 -*-

import fenics as fe

from nlhe.problems.elproblem import StrainDensityFunction

__author__ = "Alessio Nava"
__copyright__ = "2017 Alessio Nava"
__credits__ = ["Alessio Nava, the FEniCS Team"]
__license__ = "GPL"
__version__ = "3.0"
__maintainer__ = "Alessio Nava"
__email__ = "alessio.nava2@mail.polimi.it"
__status__ = "Development"


def sk_strain_density_function(ee, ey, nu):
    mu = ey / (2. * (1 + nu))
    lambda_ = ey * nu / ((1 + nu) * (1 - 2 * nu))
    return lambda_ / 2. * (fe.tr(ee)) ** 2. + mu * fe.tr(ee * ee)


strain_density_function = StrainDensityFunction(sk_strain_density_function)
