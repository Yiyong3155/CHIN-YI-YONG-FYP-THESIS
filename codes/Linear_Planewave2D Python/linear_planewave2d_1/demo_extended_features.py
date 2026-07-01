"""
Demo of Extended Linear Solver Features for FYP
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

# Setup (same as before)
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

# Material
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

# Create model
model = LinearSpectral2D(
    element._e, mesh._cpp_object, mt_facet._cpp_object,
    c0._cpp_object, rho0._cpp_object,
    sourceFrequency, sourceAmplitude, speedOfSound)

model.init()
u_n = model.u_sol()

# Create extended API wrapper
extended = LinearSpectral2DExtended(model, mesh, u_n)

if mpi_rank == 0:
    print("="*70)
    print("EXTENDED LINEAR SOLVER - FYP DEMO")
    print("="*70)

# Define regions of interest
region_near_source = {
    'xmin': 0.0, 'xmax': 0.03,
    'ymin': -0.035, 'ymax': 0.035
}

region_far_field = {
    'xmin': 0.06, 'xmax': 0.09,
    'ymin': -0.035, 'ymax': 0.035
}

# Feature 1: Write VTX with region metadata
if mpi_rank == 0:
    print("\n🎯 Feature 1: Region-aware VTX output")
extended.write_region_vtx(
    "output_with_region.bp",
    region_near_source,
    startTime, finalTime, timeStepSize, output_interval
)

# Feature 3: Generate ParaView state file
if mpi_rank == 0:
    print("\n🎯 Feature 3: Auto-generate ParaView state")
extended.save_paraview_state(
    "output_with_region.bp/data.vtx",
    region_near_source,
    "auto_clipped_view.pvsm"
)

# Feature 4: Extract region statistics
if mpi_rank == 0:
    print("\n🎯 Feature 4: Region statistics")
stats = extended.get_region_stats(region_near_source)
if mpi_rank == 0 and stats:
    print(f"Near source region stats:")
    print(f"  Max pressure: {stats['max']:.2e} Pa")
    print(f"  Min pressure: {stats['min']:.2e} Pa")
    print(f"  Mean pressure: {stats['mean']:.2e} Pa")
    print(f"  Std deviation: {stats['std']:.2e} Pa")

# Feature 2: Compare two regions
if mpi_rank == 0:
    print("\n🎯 Feature 2: Multi-region comparison")
extended.compare_regions(
    region_near_source,
    region_far_field,
    "Near Source",
    "Far Field"
)

if mpi_rank == 0:
    print("\n" + "="*70)
    print("✅ ALL FYP FEATURES DEMONSTRATED!")
    print("="*70)
    print("\n📁 Generated files:")
    print("  - output_with_region.bp/        (VTX output)")
    print("  - output_with_region.bp_region.json  (Region metadata)")
    print("  - auto_clipped_view.pvsm        (ParaView state file)")
    print("\n💡 Usage:")
    print("  paraview --state=auto_clipped_view.pvsm")
    print("="*70)

