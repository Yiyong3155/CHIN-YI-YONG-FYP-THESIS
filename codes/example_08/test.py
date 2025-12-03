from Linear import LinearSpectral2D
import basix
from dolfinx.fem import functionspace
from dolfinx.mesh import create_rectangle, CellType
from mpi4py import MPI
import numpy as np

lagrange = basix.create_element(
    basix.ElementFamily.P, 
    basix.CellType.triangle, 
    1, 
    basix.LagrangeVariant.equispaced,
    dtype=np.float32
)

# print(help(lagrange._e))

domain = create_rectangle(
    comm=MPI.COMM_WORLD,
    points=((0.0, 0.0), (2.0, 1.0)),
    n=(32, 16),
    cell_type=CellType.triangle,
    dtype=np.float32
)

LinearSpectral2D(lagrange._e, domain._cpp_object)