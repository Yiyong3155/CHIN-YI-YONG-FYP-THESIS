from Linear import LinearSpectral2D
import basix
import numpy as np

lagrange = basix.create_element(
    basix.ElementFamily.P, 
    basix.CellType.triangle, 
    1, 
    basix.LagrangeVariant.equispaced,
    dtype=np.float32
)

print(help(lagrange._e))