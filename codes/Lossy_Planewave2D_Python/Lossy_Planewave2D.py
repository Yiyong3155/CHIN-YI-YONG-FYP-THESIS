"""
Linear solver for the 2D planewave problem with attenuation
- structured mesh
- first-order Sommerfeld ABC
- submesh VTX output on custom region
===========================================================
Copyright (C) 2022 Adeeb Arif Kor
"""
 
import numpy as np
from mpi4py import MPI
import basix, basix.ufl
from dolfinx.common import Timer
from dolfinx.fem import Function, functionspace
from dolfinx.io import XDMFFile
from dolfinx import cpp
from fusx import LossySpectral2D, compute_diffusivity_of_sound
 
mpi_rank = MPI.COMM_WORLD.rank
mpi_size = MPI.COMM_WORLD.size
 
# =====================================================================
# SOURCE PARAMETERS
# =====================================================================
sourceFrequency = 0.5e6
sourceAmplitude = 60000.0
period = 1.0 / sourceFrequency
angularFrequency = 2.0 * np.pi * sourceFrequency
 
# =====================================================================
# MATERIAL PARAMETERS
# =====================================================================
speedOfSound = 1500.0
density = 1000.0
attenuationCoefficientdB = 50.0
attenuationCoefficientNp = attenuationCoefficientdB / 20.0 * np.log(10.0)
diffusivityOfSound = compute_diffusivity_of_sound(
    angularFrequency, speedOfSound, attenuationCoefficientNp
)
 
# =====================================================================
# DOMAIN PARAMETERS
# =====================================================================
domainLength = 0.12
 
# =====================================================================
# FINITE ELEMENT PARAMETERS
# =====================================================================
degreeOfBasis = 4
 
# ── Custom output region (metres) ────────────────────────────────────
REGION_MIN = [-0.06, 0.00]   # [x_min, y_min]
REGION_MAX = [ 0.06, 0.12]   # [x_max, y_max]
GRID_NX    = 151
GRID_NY    = 151
# Output: vtx_output/lossy_region.bp
# ─────────────────────────────────────────────────────────────────────
 
# =====================================================================
# READ MESH
# =====================================================================
with XDMFFile(MPI.COMM_WORLD, "mesh.xdmf", "r") as fmesh:
    mesh = fmesh.read_mesh(name="planewave_2d_1")
    mesh.topology.create_entities(1)
    mt_cell = fmesh.read_meshtags(mesh, name="planewave_2d_1_cells")
    mt_facet = fmesh.read_meshtags(mesh, name="planewave_2d_1_facets")
 
mesh.topology.create_connectivity(1, 2)
 
# =====================================================================
# MESH SIZE
# =====================================================================
tdim = mesh.topology.dim
num_cells = mesh.topology.index_map(tdim).size_local
cell_indices = np.arange(num_cells, dtype=np.int32)
h = cpp.mesh.h(mesh._cpp_object, tdim, cell_indices)
meshSizeMinLocal = np.min(h)
meshSizeMinGlobal = MPI.COMM_WORLD.allreduce(meshSizeMinLocal, op=MPI.MIN)
 
# =====================================================================
# FINITE ELEMENT SPACE
# =====================================================================
family = basix.ElementFamily.P
variant = basix.LagrangeVariant.gll_warped
cell_type = mesh.basix_cell()
element = basix.ufl.element(family, cell_type, degreeOfBasis, variant)
 
element_DG = basix.ufl.element(
    basix.ElementFamily.P, mesh.basix_cell(), 0,
    basix.LagrangeVariant.gll_warped, discontinuous=True
)
V_DG = functionspace(mesh, element_DG)
 
# =====================================================================
# MATERIAL PROPERTIES
# =====================================================================
c0     = Function(V_DG)
rho0   = Function(V_DG)
delta0 = Function(V_DG)
 
cells_1 = mt_cell.find(1)
c0.x.array[cells_1]     = speedOfSound;       c0.x.scatter_forward()
rho0.x.array[cells_1]   = density;            rho0.x.scatter_forward()
delta0.x.array[cells_1] = diffusivityOfSound; delta0.x.scatter_forward()
 
# =====================================================================
# TEMPORAL PARAMETERS
# =====================================================================
CFL            = 0.4
timeStepSize   = CFL * meshSizeMinGlobal / (speedOfSound * degreeOfBasis**2)
stepPerPeriod  = int(period / timeStepSize) + 1
timeStepSize   = period / stepPerPeriod
startTime      = 0.0
finalTime      = domainLength / speedOfSound + 4.0 / sourceFrequency
numberOfStep   = int((finalTime - startTime) / timeStepSize) + 1
output_interval = stepPerPeriod   # write once per acoustic period
 
# =====================================================================
# PRINT PARAMETERS
# =====================================================================
if mpi_rank == 0:
    print("=" * 70)
    print("PLANEWAVE 2D SIMULATION WITH ATTENUATION (Lossy) - SUBMESH OUTPUT")
    print("=" * 70)
    print(f"Speed of sound:       {speedOfSound} m/s")
    print(f"Density:              {density} kg/m³")
    print(f"Attenuation:          {attenuationCoefficientdB} dB/m")
    print(f"Diffusivity:          {diffusivityOfSound:.6e} m²/s")
    print(f"Source frequency:     {sourceFrequency/1e6} MHz")
    print(f"Source amplitude:     {sourceAmplitude/1000} kPa")
    print(f"Domain length:        {domainLength*1000} mm")
    print(f"Polynomial degree:    {degreeOfBasis}")
    print(f"Min mesh size:        {meshSizeMinGlobal*1e6:.2f} μm")
    print(f"Time step:            {timeStepSize*1e9:.3f} ns")
    print(f"Steps/period:         {stepPerPeriod}")
    print(f"Total steps:          {numberOfStep}")
    print(f"Output interval:      every {output_interval} steps (once per period)")
    print(f"Region:               x={REGION_MIN[0]} → {REGION_MAX[0]} m")
    print(f"                      y={REGION_MIN[1]} → {REGION_MAX[1]} m")
    print(f"Grid:                 {GRID_NX} x {GRID_NY}")
    print(f"Output:               vtx_output/lossy_region.bp")
    print("=" * 70)
 
# =====================================================================
# CREATE MODEL WITH CUSTOM REGION
# =====================================================================
model = LossySpectral2D(
    element.basix_element._e,
    mesh._cpp_object,
    mt_facet._cpp_object,
    c0._cpp_object,
    rho0._cpp_object,
    delta0._cpp_object,
    sourceFrequency,
    sourceAmplitude,
    speedOfSound,
    region_min=REGION_MIN,
    region_max=REGION_MAX,
    grid_nx=GRID_NX,
    grid_ny=GRID_NY,
)
 
model.init()
 
# =====================================================================
# SOLVE
# =====================================================================
if mpi_rank == 0:
    print("\n🚀 Starting lossy simulation...")
 
with Timer() as tsolve:
    model.rk4(startTime, finalTime, timeStepSize, output_interval, "lossy_region")
 
if mpi_rank == 0:
    print(f"\n{'='*70}")
    print("✅ LOSSY SIMULATION COMPLETE!")
    print(f"{'='*70}")
    print(f"Solve time:   {tsolve.elapsed()[0]:.2f} s")
    print(f"DOFs:         {model.number_of_dofs()}")
    print(f"Output:       vtx_output/lossy_region.bp")
    print(f"ParaView:     File > Open > vtx_output/lossy_region.bp")
    print(f"{'='*70}")
 