"""
Linear solver for the 2D circular penetrable scatterer
- circular scatterer
- wavelength < scatterer radius
- submesh VTX output on custom region (animation ready)
======================================================
Copyright (C) 2024 Adeeb Arif Kor
"""
 
import numpy as np
from mpi4py import MPI
import basix, basix.ufl
from dolfinx.common import Timer
from dolfinx.fem import Function, functionspace
from dolfinx.io import XDMFFile
from dolfinx import cpp
from fusx import LinearSpectral2D
 
mpi_rank = MPI.COMM_WORLD.rank
mpi_size = MPI.COMM_WORLD.size
 
# =====================================================================
# MATERIAL PARAMETERS
# =====================================================================
speedOfSound1 = 1500.0   # m/s (medium 1)
density1      = 1000.0   # kg/m³ (medium 1)
speedOfSound2 = 3500.0   # m/s (scatterer)
density2      = 1900.0   # kg/m³ (scatterer)
 
# =====================================================================
# SOURCE PARAMETERS
# =====================================================================
sourceFrequency = 2000.0
sourceSpeed     = 1.0
sourceAmplitude = density1 * speedOfSound1 * sourceSpeed
period          = 1.0 / sourceFrequency
 
# =====================================================================
# DOMAIN PARAMETERS
# =====================================================================
wavelength      = speedOfSound1 / sourceFrequency
scattererRadius = 3.0 * wavelength
domainScale     = 10.0
domainLength    = 1.5 * domainScale * scattererRadius + 2 * scattererRadius
 
# =====================================================================
# FINITE ELEMENT PARAMETERS
# =====================================================================
degreeOfBasis = 4
 
# ── Custom output region (metres) ────────────────────────────────────
# Adjust these to zoom into the area of interest around the scatterer
REGION_MIN = [-scattererRadius * 2, -scattererRadius * 2]
REGION_MAX = [ scattererRadius * 2,  scattererRadius * 2]
GRID_NX    = 251
GRID_NY    = 251
# Output: vtx_output/scatterer_region.bp
# ─────────────────────────────────────────────────────────────────────
 
# =====================================================================
# READ MESH
# =====================================================================
with XDMFFile(MPI.COMM_WORLD, "../mesh.xdmf", "r") as fmesh:
    mesh = fmesh.read_mesh(name="circular_scatterer_penetrable_2d_1")
    mesh.topology.create_entities(1)
    mt_cell  = fmesh.read_meshtags(mesh, name="circular_scatterer_penetrable_2d_1_cells")
    mt_facet = fmesh.read_meshtags(mesh, name="circular_scatterer_penetrable_2d_1_facets")
 
mesh.topology.create_connectivity(1, 2)
 
# =====================================================================
# MESH SIZE
# =====================================================================
tdim        = mesh.topology.dim
num_cells   = mesh.topology.index_map(tdim).size_local
cell_indices = np.arange(num_cells, dtype=np.int32)
h = cpp.mesh.h(mesh._cpp_object, tdim, cell_indices)
meshSizeMinLocal  = np.min(h)
meshSizeMinGlobal = MPI.COMM_WORLD.allreduce(meshSizeMinLocal, op=MPI.MIN)
 
# =====================================================================
# FINITE ELEMENT SPACE
# =====================================================================
element = basix.ufl.element(
    basix.ElementFamily.P, mesh.basix_cell(), degreeOfBasis,
    basix.LagrangeVariant.gll_warped
)
element_DG = basix.ufl.element(
    basix.ElementFamily.P, mesh.basix_cell(), 0,
    basix.LagrangeVariant.gll_warped, discontinuous=True
)
V_DG = functionspace(mesh, element_DG)
 
# =====================================================================
# MATERIAL PROPERTIES
# =====================================================================
c0   = Function(V_DG)
rho0 = Function(V_DG)
 
cells_1 = mt_cell.find(1)   # medium
cells_2 = mt_cell.find(2)   # scatterer
 
c0.x.array[cells_1] = speedOfSound1
c0.x.array[cells_2] = speedOfSound2
c0.x.scatter_forward()
 
rho0.x.array[cells_1] = density1
rho0.x.array[cells_2] = density2
rho0.x.scatter_forward()
 
# =====================================================================
# TEMPORAL PARAMETERS
# =====================================================================
CFL            = 0.9
timeStepSize   = CFL * meshSizeMinGlobal / (speedOfSound2 * degreeOfBasis**2)
stepPerPeriod  = int(period / timeStepSize) + 1
timeStepSize   = period / stepPerPeriod
startTime      = 0.0
finalTime      = domainLength / speedOfSound1 + 8.0 / sourceFrequency
numberOfStep   = int((finalTime - startTime) / timeStepSize) + 1
output_interval = stepPerPeriod   # one frame per acoustic period
 
# =====================================================================
# PRINT PARAMETERS
# =====================================================================
if mpi_rank == 0:
    print("=" * 70)
    print("CIRCULAR SCATTERER 2D (Linear) - SUBMESH VTX OUTPUT")
    print("=" * 70)
    print(f"Scatterer radius:          {scattererRadius:.4f} m")
    print(f"Wavelength:                {wavelength:.4f} m")
    print(f"Speed of sound (medium):   {speedOfSound1} m/s")
    print(f"Speed of sound (scatter):  {speedOfSound2} m/s")
    print(f"Density (medium):          {density1} kg/m³")
    print(f"Density (scatter):         {density2} kg/m³")
    print(f"Source frequency:          {sourceFrequency} Hz")
    print(f"Source amplitude:          {sourceAmplitude:.2f} Pa")
    print(f"Domain length:             {domainLength:.4f} m")
    print(f"Polynomial degree:         {degreeOfBasis}")
    print(f"Min mesh size:             {meshSizeMinGlobal:.2e} m")
    print(f"Time step:                 {timeStepSize:.6e} s")
    print(f"Steps/period:              {stepPerPeriod}")
    print(f"Total steps:               {numberOfStep}")
    print(f"Output interval:           every {output_interval} steps (once per period)")
    print(f"Region:                    x={REGION_MIN[0]:.3f} → {REGION_MAX[0]:.3f} m")
    print(f"                           y={REGION_MIN[1]:.3f} → {REGION_MAX[1]:.3f} m")
    print(f"Grid:                      {GRID_NX} x {GRID_NY}")
    print(f"Output:                    vtx_output/scatterer_region.bp")
    print("=" * 70)
 
# =====================================================================
# CREATE MODEL WITH CUSTOM REGION
# =====================================================================
model = LinearSpectral2D(
    element.basix_element._e,
    mesh._cpp_object,
    mt_facet._cpp_object,
    c0._cpp_object,
    rho0._cpp_object,
    sourceFrequency,
    sourceAmplitude,
    speedOfSound1,
    region_min=REGION_MIN,
    region_max=REGION_MAX,
    grid_nx=GRID_NX,
    grid_ny=GRID_NY,
)
 
model.init()
 
# =====================================================================
# SOLVE WITH ANIMATION OUTPUT
# =====================================================================
if mpi_rank == 0:
    print("\n🚀 Starting simulation...")
 
with Timer() as tsolve:
    model.rk4(startTime, finalTime, timeStepSize, output_interval, "scatterer_region")
 
if mpi_rank == 0:
    print(f"\n{'='*70}")
    print("✅ SIMULATION COMPLETE!")
    print(f"{'='*70}")
    print(f"Solve time:  {tsolve.elapsed()[0]:.2f} s")
    print(f"Time/step:   {tsolve.elapsed()[0]/numberOfStep:.4f} s")
    print(f"DOFs:        {model.number_of_dofs()}")
    print(f"Output:      vtx_output/scatterer_region.bp")
    print(f"ParaView:    File > Open > vtx_output/scatterer_region.bp")
    print(f"             → colour by 'u' → press Play for animation")
    print(f"{'='*70}")
 