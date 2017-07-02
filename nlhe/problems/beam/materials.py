#!/usr/bin/env python
# -*- coding: utf-8 -*-

from nlhe.problems.elproblem import Materials

__author__ = "Alessio Nava"
__copyright__ = "2017 Alessio Nava"
__credits__ = ["Alessio Nava, the FEniCS Team"]
__license__ = "GPL"
__version__ = "3.0"
__maintainer__ = "Alessio Nava"
__email__ = "alessio.nava2@mail.polimi.it"
__status__ = "Development"

materials = Materials()

materials.add_material(1, 210e9, 0.33)
