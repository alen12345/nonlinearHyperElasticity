#!/usr/bin/env python
# -*- coding: utf-8 -*-

import importlib
import os
import errno

from nlhe.solvers import nlhe

__author__ = "Alessio Nava"
__copyright__ = "2017 Alessio Nava"
__credits__ = ["Alessio Nava, the FEniCS Team"]
__license__ = "GPL"
__version__ = "3.0"
__maintainer__ = "Alessio Nava"
__email__ = "alessio.nava2@mail.polimi.it"
__status__ = "Development"


class Materials:
    def __init__(self):
        self._materials = {}

    def get_materials(self):
        return self._materials

    def add_material(self, mid, nu, eyoung):
        self._materials[mid] = [nu, eyoung, nu]


class BoundaryConditions:
    def __init__(self):
        self._bondary_conditions = {}

    def get_boundary_conditions(self):
        return self._bondary_conditions

    def add_boundary_condition(self, bid, bc):
        self._bondary_conditions[bid] = bc


class SurfaceLoads:
    def __init__(self):
        self._surface_loads = {}

    def get_surface_loads(self):
        return self._surface_loads

    def add_surface_load(self, sid, sload):
        self._surface_loads[sid] = sload


class BodyLoads:
    def __init__(self):
        self._body_loads = {}

    def get_body_loads(self):
        return self._body_loads

    def add_body_loads(self, bid, bload):
        self._body_loads[bid] = bload


class StrainDensityFunction:
    def __init__(self, strain_density_function):
        self._strain_density_function = strain_density_function

    def get_strain_density_function(self):
        return self._strain_density_function


class ElasticProblem:
    """
    Provides a consistent interface for the definition of a solid mechanics problem. The interface
    makes no assumption on the type of problem or solver implemented; the only mandatory requirement
    is to supply a valid **strain energy density function**.
    """

    def _path_setup(self):
        """
        First function called by the constructor to define the relevant paths to define
        the elastic problem.

        :return: no value is returned, only local private variables are affected.
        :rtype: None
        """
        self._mesh_path = os.path.join('nlhe', 'gmsh', self._analysis_name, self._analysis_name + '.geo')
        self._dolfin_mesh_dir = os.path.join('nlhe', 'meshes', self._analysis_name)
        self._problem_dir = os.path.join('nlhe', 'problems', self._analysis_name)
        self._materials_path = os.path.join(self._problem_dir, 'materials.py')
        self._bcs_path = os.path.join(self._problem_dir, 'bcs.py')
        self._sloads_path = os.path.join(self._problem_dir, 'sloads.py')
        self._bloads_path = os.path.join(self._problem_dir, 'bloads.py')
        self._sdf_path = os.path.join(self._problem_dir, 'sdf.py')

    def _mesh_setup(self):
        """
        Checks if the analysis that the user wants to run is in the specified path. If the path to
        the converted dolfin mesh has not been created yet, a new directory is added.

        :return: does not return any values.
        :rtype: None
        """
        if not os.path.isfile(self._mesh_path):
            raise IOError(errno.ENOENT, os.strerror(errno.ENOENT), self._mesh_path)
        if not os.path.isdir(self._dolfin_mesh_dir):
            os.makedirs(self._dolfin_mesh_dir)

    def _problem_setup(self):
        """
        Checks if all the ingredients of the problem have been defined by the user in the correct path.
        Otherwise, an *exception* indicating the right path is raised.

        :return: does not return any values.
        :rtype: None
        """
        if not os.path.isdir(self._problem_dir):
            raise IOError(errno.ENOENT, "Project structure for the problem is not respected. Check the docs.")
        if not os.path.isfile(self._materials_path):
            raise IOError(errno.ENOENT, "Cannot find the specifications for the materials" +
                          ", which should be located at " + self._materials_path)
        if not os.path.isfile(self._bcs_path):
            raise IOError(errno.ENOENT, "Cannot find the specifications for the boundary conditions" +
                          ", which should be located at " + self._bcs_path)
        if not os.path.isfile(self._sloads_path):
            raise IOError(errno.ENOENT, "Cannot find the specifications for the surface loads" +
                          ", which should be located at " + self._sloads_path)
        if not os.path.isfile(self._bloads_path):
            raise IOError(errno.ENOENT, "Cannot find the specifications for the body loads" +
                          ", which should be located at " + self._bloads_path)
        if not os.path.isfile(self._sdf_path):
            raise IOError(errno.ENOENT, "Cannot find the specifications for the strain density function" +
                          ", which should be located at " + self._sdf_path)

    def _material_setup(self):
        """
        Dynamically imports the user supplied module containing the definition of the materials.

        :return: does not return any values.
        :rtype: None
        """
        materials_module = importlib.import_module('nlhe.problems.' + self._analysis_name + '.materials')
        self._materials = materials_module.materials.get_materials()

    def _boundary_conditions_setup(self):
        """
        Dynamically imports the user supplied module containing the definition of the boundary conditions.

        :return: does not return any values.
        :rtype: None
        """
        bcs_module = importlib.import_module('nlhe.problems.' + self._analysis_name + '.bcs')
        self._boundary_conditions = bcs_module.boundary_conditions.get_boundary_conditions()

    def _surface_loads_setup(self):
        """
        Dynamically imports the user supplied module containing the definition of the surface loads.

        :return: does not return any values.
        :rtype: None
        """
        surface_loads_module = importlib.import_module('nlhe.problems.' + self._analysis_name + '.sloads')
        self._surface_loads = surface_loads_module.surface_loads.get_surface_loads()

    def _body_loads_setup(self):
        """
        Dynamically imports the user supplied module containing the definition of the body loads.

        :return: does not return any values.
        :rtype: None
        """
        body_loads_module = importlib.import_module('nlhe.problems.' + self._analysis_name + '.bloads')
        self._body_loads = body_loads_module.body_loads.get_body_loads()

    def _strain_density_function_setup(self):
        """
        Dynamically imports the user supplied module containing the definition of the
        strain energy density function.

        :return: does not return any values.
        :rtype: None
        """
        sdf_module = importlib.import_module('nlhe.problems.' + self._analysis_name + '.sdf')
        self._strain_density_function = sdf_module.strain_density_function.get_strain_density_function()

    def solve(self):
        """
        Invokes the default solver, which for the moment is a general purpose isotropic nonlinear
        hyperelasticity solver.

        :return: does not return any values.
        :rtype: None
        :note: For the time being, there is only one default solver implemented.
        """
        nlhe.nlhe_solver(self._analysis_name, self._materials, self._boundary_conditions,
                         self._surface_loads, self._body_loads,
                         self._strain_density_function)

    def __init__(self, analysis_name):
        """
        The default constructor sets up the problem as described in the auxiliary function.
        It is up to the user to *explicitly invoke the solver*.

        :param analysis_name: name of the analysis passed on the commandline and parsed.
        :type analysis_name: str
        """
        self._analysis_name = analysis_name
        self._path_setup()
        self._mesh_setup()
        self._problem_setup()
        self._material_setup()
        self._boundary_conditions_setup()
        self._surface_loads_setup()
        self._body_loads_setup()
        self._strain_density_function_setup()
