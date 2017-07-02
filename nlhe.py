#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

import argparse

from nlhe.problems import elproblem

__author__ = "Alessio Nava"
__copyright__ = "2017 Alessio Nava"
__credits__ = ["Alessio Nava, the FEniCS Team"]
__license__ = "GPL"
__version__ = "3.0"
__maintainer__ = "Alessio Nava"
__email__ = "alessio.nava2@mail.polimi.it"
__status__ = "Development"


def main(args):
    elastic_problem = elproblem.ElasticProblem(args.analysis_name)
    elastic_problem.solve()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Solve the nonlinear hyperelastic equation.',
                                     epilog="""Pre and post-procesing capabilties with gmsh and paraview, solver in 
                                     pure Python thanks to FEniCS. """
                                     )
    parser.add_argument('analysis_name', type=str,
                        help="""name of the analysis to be executed.""")
    parser.add_argument('-s', '--solver', default='petsc', type=str,
                        help="""See the code for the available solvers.""")
    parser.add_argument('--version', action='version', version='%(prog)s 1.0.0')

    arguments = parser.parse_args()
    main(arguments)
