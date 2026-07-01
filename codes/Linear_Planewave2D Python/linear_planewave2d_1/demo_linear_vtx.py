"""
Linear solver with submesh VTX output on custom region
Output stored in: ./vtx_output/pressure_region.bp
"""
import numpy as np
from mpi4py import MPI
import basix
import basix.ufl
from dolfinx.common import Timer
from dolfinx.fem import Function, functionspace
from dolfinx.io import XDMFFile
from dolfinx import cpp
from fusx import LinearSpectral2D

mpi_rank = MPI.COMM_WORLD.rank

# Parameters
sourceFrequency = 0.5e6
sourceAmplitude  = 60000
speedOfSound     = 1500
density          = 1000
domainLength     = 0.12
degreeOfBasis    = 4

# ── Custom output region (metres) ──────────────────────────────────────────
REGION_MIN = [-0.06, 0.00]   # [x_min, y_min]
REGION_MAX = [ 0.06, 0.12]   # [x_max, y_max]
GRID_NX    = 151
GRID_NY    = 151
# Output: ./vtx_output/pressure_region.bp
# ───────────────────────────────────────────────────────────────────────────

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

V_DG  = functionspace(mesh, ("DG", 0))
c0    = Function(V_DG); c0.x.array[:] = speedOfSound
rho0  = Function(V_DG); rho0.x.array[:] = density

period          = 1 / sourceFrequency
CFL             = 0.9
timeStepSize    = CFL * meshSize / (speedOfSound * degreeOfBasis**2)
stepPerPeriod   = int(period / timeStepSize + 1)
timeStepSize    = period / stepPerPeriod
startTime       = 0.0
finalTime       = domainLength / speedOfSound + 4.0 / sourceFrequency
output_interval = stepPerPeriod

if mpi_rank == 0:
    print("=" * 70)
    print("LINEAR SOLVER - SUBMESH VTX OUTPUT")
    print("=" * 70)
    print(f"Region:      x=[{REGION_MIN[0]}, {REGION_MAX[0]}]  y=[{REGION_MIN[1]}, {REGION_MAX[1]}]  (metres)")
    print(f"Grid:        {GRID_NX} x {GRID_NY} = {GRID_NX*GRID_NY} points")
    print(f"Output:      vtx_output/pressure_region.bp")
    print(f"Interval:    every {output_interval} steps (once per period)")
    print(f"Final time:  {finalTime*1e6:.2f} µs  |  dt: {float(timeStepSize)*1e9:.4f} ns")
    print("=" * 70)

model = LinearSpectral2D(
    element._e,
    mesh._cpp_object,
    mt_facet._cpp_object,
    c0._cpp_object,
    rho0._cpp_object,
    sourceFrequency,
    sourceAmplitude,
    speedOfSound,
    region_min=REGION_MIN,
    region_max=REGION_MAX,
    grid_nx=GRID_NX,
    grid_ny=GRID_NY,
)

model.init()

if mpi_rank == 0:
    print("\n🚀 Starting simulation...")

with Timer() as tsolve:
    model.rk4(startTime, finalTime, timeStepSize, output_interval, "pressure_region")

if mpi_rank == 0:
    print(f"\n{'='*70}")
    print("✅ DONE")
    print(f"Solve time:  {tsolve.elapsed()[0]:.2f} s")
    import os
    print(f"Output:      /home/shared/example_08/vtx_output/pressure_region.bp")
    print(f"ParaView:    File > Open > vtx_output/pressure_region.bp")
    print(f"{'='*70}")
