#!/usr/bin/env python
# -*- coding: utf-8 -*-

# MODULE IMPORT
from __future__ import print_function
import subprocess
from termcolor import cprint

import fenics as fe

__author__ = "Alessio Nava"
__copyright__ = "2017 Alessio Nava"
__credits__ = ["Alessio Nava, the FEniCS Team"]
__license__ = "GPL"
__version__ = "3.0"
__maintainer__ = "Alessio Nava"
__email__ = "alessio.nava2@mail.polimi.it"
__status__ = "Development"


def nlhe_solver(mesh_name, materials, boundary_conditions,
                surface_loads, body_forces,
                strain_density_function):
    """
    **Nonlinear hyperelasticity solver**, with pre and post capabiltities.

    The problem is treated in large deformation, proposing a generic hyperelastic model, with
    isotropic beahaviour. Pre and post capabilities are provided through *gmsh and paraview*.

    .. math::
        :label: integral

        \int f(x) \mathrm{d}x

    :param mesh_name: The mesh name automatically retrieves the project structure.
    :param materials: Dictionary of materials.
    :param boundary_conditions: Dictionary o boundary conditions.
    :param surface_loads: Dictionary of surface loads.
    :param body_forces: Dictionary
    :param strain_density_function: Function describing the material model.
    :type mesh_name: string
    :return: Does not return a value, but writes the results for further postprocessing.
    :rtype: None
    :warning: No explicit warning here.
    :note: This is crazy, isn't it?
    :seealso: What should I also look at?
    """
    # COMPILER OPTIONS
    fe.parameters.linear_algebra_backend = 'PETSc'
    fe.parameters.form_compiler.cpp_optimize_flags = '-O3'
    fe.parameters.form_compiler.quadrature_degree = 3

    # IMPORT MESH
    cprint("\nMESH PRE-PROCESSING", 'blue', attrs=['bold'])
    # Import mesh and groups
    cprint("Creating gmsh mesh...", 'green')
    subprocess.check_output("gmsh nlhe/gmsh" + 2 * ("/" + mesh_name) + ".geo -3", shell=True)
    cprint("Converting mesh to DOLFIN format...", 'green')
    subprocess.check_output("dolfin-convert nlhe/gmsh" + 2 * ("/" + mesh_name) +
                            ".msh nlhe/meshes" + 2 * ("/" + mesh_name) + ".xml", shell=True)
    cprint("Importing mesh in FEniCS...", 'green')
    mesh = fe.Mesh("nlhe/meshes" + 2 * ("/" + mesh_name) + ".xml")
    cprint("Generating boundaries and subdomains...", 'green')
    subdomains = fe.MeshFunction("size_t", mesh,
                                 "nlhe/meshes" + 2 * ("/" + mesh_name) + "_physical_region.xml")
    boundaries = fe.MeshFunction("size_t", mesh,
                                 "nlhe/meshes" + 2 * ("/" + mesh_name) + "_facet_region.xml")

    # Redefine the integration measures
    dxp = fe.Measure('dx', domain=mesh, subdomain_data=subdomains)
    dsp = fe.Measure('ds', domain=mesh, subdomain_data=boundaries)

    # FINITE ELEMENT SPACES
    # Finite element spaces
    W = fe.FunctionSpace(mesh, 'P', 1)
    V = fe.VectorFunctionSpace(mesh, 'P', 1)
    Z = fe.TensorFunctionSpace(mesh, 'P', 1)

    # Finite element functions
    du = fe.TrialFunction(V)
    v = fe.TestFunction(V)
    u = fe.Function(V)

    # BOUNDARY CONDITIONS
    # Collect Dirichlet conditions
    bcs = []
    for i in boundary_conditions:
        bc = fe.DirichletBC(V, boundary_conditions[i],
                            boundaries, i)
    bcs.append(bc)

    # LOAD DEFINITIONS
    # Collect the surface loads
    integrals_S = []
    for i in surface_loads:
        gLoad = surface_loads[i]
        integrals_S.append(fe.dot(gLoad, u) * dsp(i))

    # Collect the volume loads
    integrals_V = []
    for i in body_forces:
        fLoad = body_forces[i]
        integrals_V.append(fe.dot(fLoad, u) * dxp(i))

    # VARIATIONAL FORMULATION
    # Large displacements kinematics
    def large_kinematics(u):
        d = u.geometric_dimension()
        I = fe.Identity(d)
        Fgrad = I + fe.grad(u)
        C = Fgrad.T * Fgrad
        E = 0.5 * (C - I)
        return E

    E = large_kinematics(u)
    E = fe.variable(E)

    # Collect the strain density functions
    integrals_E = []
    for i in materials:
        psi = strain_density_function(E, materials[i][0], materials[i][1])
        integrals_E.append(psi * dxp(i))

    # Total potential energy
    Pi = sum(integrals_E) - sum(integrals_V) - sum(integrals_S)

    # Compute 1st variation of Pi (directional derivative about u in dir. of v)
    F = fe.derivative(Pi, u, v)

    # Compute Jacobian of F
    J = fe.derivative(F, u, du)

    # SOLVER PARAMS
    # Define the solver params
    problem = fe.NonlinearVariationalProblem(F, u, bcs, J)
    solver = fe.NonlinearVariationalSolver(problem)
    prm = solver.parameters
    prm['nonlinear_solver'] = 'newton'
    prm['newton_solver']['linear_solver'] = 'petsc'

    prm['newton_solver']['error_on_nonconvergence'] = True
    prm['newton_solver']['absolute_tolerance'] = 1E-9
    prm['newton_solver']['relative_tolerance'] = 1E-8
    prm['newton_solver']['maximum_iterations'] = 25
    prm['newton_solver']['relaxation_parameter'] = 1.0

    prm['newton_solver']['lu_solver']['report'] = True
    prm['newton_solver']['lu_solver']['reuse_factorization'] = False
    prm['newton_solver']['lu_solver']['same_nonzero_pattern'] = False
    prm['newton_solver']['lu_solver']['symmetric'] = False

    prm['newton_solver']['krylov_solver']['error_on_nonconvergence'] = True
    prm['newton_solver']['krylov_solver']['absolute_tolerance'] = 1E-7
    prm['newton_solver']['krylov_solver']['relative_tolerance'] = 1E-5
    prm['newton_solver']['krylov_solver']['maximum_iterations'] = 1000
    prm['newton_solver']['krylov_solver']['nonzero_initial_guess'] = True

    # Invoke the solver
    cprint("\nSOLUTION OF THE NONLINEAR PROBLEM", 'blue', attrs=['bold'])
    cprint("The solution of the nonlinear system is in progress...", 'red')
    solver.solve()

    # POST-PROCESSING
    cprint("\nSOLUTION POST-PROCESSING", 'blue', attrs=['bold'])
    # Save solution to file in VTK format
    cprint("Saving displacement solution to file...", 'green')
    uViewer = fe.File("nlhe/paraview" + 2 * ("/" + mesh_name) + "_displacement.pvd")
    uViewer << u

    # Maximum and minimum displacement
    u_magnitude = fe.sqrt(fe.dot(u, u))
    u_magnitude = fe.project(u_magnitude, W)
    print('Min/Max displacement:',
          u_magnitude.vector().array().min(),
          u_magnitude.vector().array().max())

    # Computation of the large deformation strains
    cprint("Computing the deformation tensor and saving to file...", 'green')
    epsilon_u = large_kinematics(u)
    epsilon_u_project = fe.project(epsilon_u, Z)
    epsilon_viewer = fe.File("nlhe/paraview" + 2 * ("/" + mesh_name) + "_strain.pvd")
    epsilon_viewer << epsilon_u_project

    # Computation of the stresses
    cprint("Stress derivation and saving to file...", 'green')
    S = fe.diff(psi, E)
    S_project = fe.project(S, Z)
    sigma_viewer = fe.File("nlhe/paraview" + 2 * ("/" + mesh_name) + "_stress.pvd")
    sigma_viewer << S_project

    # Computation of an equivalent stress
    s = S - (1. / 3) * fe.tr(S) * fe.Identity(u.geometric_dimension())
    von_Mises = fe.sqrt(3. / 2 * fe.inner(s, s))
    von_Mises_project = fe.project(von_Mises, W)
    mises_viewer = fe.File("nlhe/paraview" + 2 * ("/" + mesh_name) + "_mises.pvd")
    mises_viewer << von_Mises_project
    print("Maximum equivalent stress:",
          von_Mises_project.vector().array().max())
