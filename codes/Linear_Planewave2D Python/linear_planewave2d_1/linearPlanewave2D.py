"""
Test efficient VTX output for Linear solver
"""
import numpy as np
from mpi4py import MPI
import basix
import basix.ufl
from dolfinx.common import Timer
from dolfinx.fem import Function, functionspace
from dolfinx.io import XDMFFile, VTXWriter
from dolfinx import cpp
from fusx import LinearSpectral2D

mpi_rank = MPI.COMM_WORLD.rank

# Parameters
sourceFrequency = 0.5e6
sourceAmplitude = 60000
speedOfSound = 1500
density = 1000
domainLength = 0.12
degreeOfBasis = 4

# Read mesh
with XDMFFile(MPI.COMM_WORLD, "mesh.xdmf", "r") as fmesh:
    mesh_name = "planewave_2d_1"
    mesh = fmesh.read_mesh(name=f"{mesh_name}")
    tdim = mesh.topology.dim
    mt_cell = fmesh.read_meshtags(mesh, name=f"{mesh_name}_cells")
    mesh.topology.create_connectivity(tdim - 1, tdim)
    mt_facet = fmesh.read_meshtags(mesh, name=f"{mesh_name}_facets")

# Mesh size
numCell = mesh.topology.index_map(tdim).size_local
hmin = np.array([cpp.mesh.h(mesh._cpp_object, tdim, np.arange(numCell)).min()])
meshSize = np.zeros(1)
MPI.COMM_WORLD.Reduce(hmin, meshSize, op=MPI.MIN, root=0)
MPI.COMM_WORLD.Bcast(meshSize, root=0)

# Element
cell_type = mesh.ufl_cell().cellname()
element = basix.ufl.element(
    basix.ElementFamily.P, cell_type, 4, basix.LagrangeVariant.gll_warped,
    dtype=np.float64
).basix_element

# Material functions
V_DG = functionspace(mesh, ("DG", 0))
c0 = Function(V_DG)
c0.x.array[:] = speedOfSound
rho0 = Function(V_DG)
rho0.x.array[:] = density

# Time parameters
period = 1 / sourceFrequency
CFL = 0.9
timeStepSize = CFL * meshSize / (speedOfSound * degreeOfBasis**2)
stepPerPeriod = int(period / timeStepSize + 1)
timeStepSize = period / stepPerPeriod
startTime = 0.0
finalTime = domainLength / speedOfSound + 4.0 / sourceFrequency

output_interval = stepPerPeriod

if mpi_rank == 0:
    print("=" * 70)
    print("EFFICIENT VTX OUTPUT TEST - LINEAR SOLVER")
    print("=" * 70)
    print(f"Output interval: {output_interval} steps")
    print("=" * 70)

# Create model
model = LinearSpectral2D(
    element._e, 
    mesh._cpp_object, 
    mt_facet._cpp_object, 
    c0._cpp_object,
    rho0._cpp_object,
    sourceFrequency,
    sourceAmplitude,
    speedOfSound)

model.init()
u_n = model.u_sol()

# Run with VTX output
if mpi_rank == 0:
    print("\n🚀 Starting simulation with VTX output...")

    with Timer() as tsolve:
        model.rk4(startTime, finalTime, timeStepSize, output_interval, "pressure_region")

if mpi_rank == 0:
    print(f"\n{'='*70}")
    print("✅ SUCCESS!")
    print(f"{'='*70}")
    print(f"Solve time: {tsolve.elapsed()[0]:.2f} seconds")
    print(f"Output: output.bp/")
    print(f"Method: Efficient C++ VTX output!")
    print(f"{'='*70}")
