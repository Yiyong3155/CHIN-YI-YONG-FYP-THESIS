import numpy as np
from mpi4py import MPI

import basix
import basix.ufl
from dolfinx.fem import Function, functionspace
from dolfinx.io import XDMFFile, VTXWriter
from dolfinx import cpp

from Linear import LinearSpectral2D

with XDMFFile(MPI.COMM_WORLD, "../mesh.xdmf", "r") as fmesh:
    mesh_name = "planewave_2d_1"
    mesh = fmesh.read_mesh(name=f"{mesh_name}")
    tdim = mesh.topology.dim
    mt_cell = fmesh.read_meshtags(mesh, name=f"{mesh_name}_cells")
    mesh.topology.create_connectivity(tdim - 1, tdim)
    mt_facet = fmesh.read_meshtags(mesh, name=f"{mesh_name}_facets")

# Define cell, finite element and function space
cell_type = mesh.ufl_cell().cellname()
element = basix.ufl.element(
    basix.ElementFamily.P, cell_type, 4, basix.LagrangeVariant.gll_warped,
    dtype=np.float64
).basix_element

# Define a DG function space for the physical parameters of the domain
V_DG = functionspace(mesh, ("DG", 0))
c0 = Function(V_DG)
c0.x.array[:] = 1.0

rho0 = Function(V_DG)
rho0.x.array[:] = 2.0

LinearSpectral2D(
    element._e, 
    mesh._cpp_object, 
    mt_facet._cpp_object, 
    c0._cpp_object,
    rho0._cpp_object,
    10.0,
    15.0,
    20.0)