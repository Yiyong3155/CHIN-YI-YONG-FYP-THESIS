"""
Demo: Export custom region in PNG and CSV formats
"""
import numpy as np
from mpi4py import MPI
import basix
import basix.ufl
from dolfinx.fem import Function, functionspace
from dolfinx.io import XDMFFile
from dolfinx import cpp
from fusx import LinearSpectral2D
from linear_extended import LinearSpectral2DExtended

mpi_rank = MPI.COMM_WORLD.rank

# Quick setup
sourceFrequency = 0.5e6
sourceAmplitude = 60000
speedOfSound = 1500
density = 1000

with XDMFFile(MPI.COMM_WORLD, "mesh.xdmf", "r") as fmesh:
    mesh_name = "planewave_2d_1"
    mesh = fmesh.read_mesh(name=f"{mesh_name}")
    tdim = mesh.topology.dim
    mt_cell = fmesh.read_meshtags(mesh, name=f"{mesh_name}_cells")
    mesh.topology.create_connectivity(tdim - 1, tdim)
    mt_facet = fmesh.read_meshtags(mesh, name=f"{mesh_name}_facets")

numCell = mesh.topology.index_map(tdim).size_local
hmin = np.array([cpp.mesh.h(mesh._cpp_object, tdim, np.arange(numCell)).min()])
meshSize = np.zeros(1)
MPI.COMM_WORLD.Reduce(hmin, meshSize, op=MPI.MIN, root=0)
MPI.COMM_WORLD.Bcast(meshSize, root=0)

cell_type = mesh.ufl_cell().cellname()
element = basix.ufl.element(
    basix.ElementFamily.P, cell_type, 4, basix.LagrangeVariant.gll_warped,
    dtype=np.float64
).basix_element

V_DG = functionspace(mesh, ("DG", 0))
c0 = Function(V_DG)
c0.x.array[:] = speedOfSound
rho0 = Function(V_DG)
rho0.x.array[:] = density

period = 1 / sourceFrequency
CFL = 0.9
timeStepSize = CFL * meshSize / (speedOfSound * 4**2)
stepPerPeriod = int(period / timeStepSize + 1)
timeStepSize = period / stepPerPeriod
startTime = 0.0
finalTime = 0.12 / speedOfSound + 4.0 / sourceFrequency

# Create model
model = LinearSpectral2D(
    element._e, mesh._cpp_object, mt_facet._cpp_object,
    c0._cpp_object, rho0._cpp_object,
    sourceFrequency, sourceAmplitude, speedOfSound)

model.init()
u_n = model.u_sol()

# Run simulation
if mpi_rank == 0:
    print("="*70)
    print("CUSTOM REGION EXPORT - PNG & CSV")
    print("="*70)
    print("\n🚀 Running simulation...")

from dolfinx.io import VTXWriter
with VTXWriter(mesh.comm, "temp_output.bp", [u_n], "bp5") as vtx:
    vtx_cpp = vtx._cpp_object
    model.rk4(startTime, finalTime, timeStepSize, stepPerPeriod, vtx_cpp)

# Create extended API
extended = LinearSpectral2DExtended(model, mesh, u_n)

# Define custom region
region = {
    'xmin': 0.0, 'xmax': 0.03,
    'ymin': -0.035, 'ymax': 0.035
}

# Export as CSV (for ParaView or analysis)
extended.export_region_csv(region, "custom_region.csv")

if mpi_rank == 0:
    print("\n" + "="*70)
    print("✅ EXPORT COMPLETE!")
    print("="*70)
    print("\n📁 Output files:")
    print("  - custom_region.csv        (Custom region data)")
    print("\n💡 Use in ParaView:")
    print("  1. Open custom_region.csv")
    print("  2. Filters → Table To Points")
    print("  3. X Column: x, Y Column: y, Z Column: 0")
    print("  4. Color by: pressure")
    print("\n💡 Use in Python/Matlab:")
    print("  import pandas as pd")
    print("  df = pd.read_csv('custom_region.csv')")
    print("  # Now analyze df['pressure']")
    print("="*70)

