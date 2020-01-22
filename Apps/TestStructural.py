from GridLib import *
from FieldsLib import *
from LinearSystemLib import *
from GeoLib import *
import numpy as np
import pylab as pl

# -------------- GRID DATA ----------------------------
L_0 = 4.
L_1 = 6.
L = L_0 + L_1
nVertices = 15
nodesCoord, elemConn = createGridData(L, nVertices)
gridData = GridData()
gridData.setElementConnectivity(elemConn)
gridData.setNodeCoordinates(nodesCoord)
grid = Grid_1D(gridData)
# -----------------------------------------------------

# -------------- PROPERTIES ----------------------------
M = 1.3e8
modulus = ScalarField(grid.getNumberOfRegions())
modulus.setValue(grid.getRegions()[0], M)
rho = 2300.
density = ScalarField(grid.getNumberOfRegions())
density.setValue(grid.getRegions()[0], rho)
g = -9.81
# -----------------------------------------------------

# -------------- NUMERICAL SOLUTION -------------------
ls = LinearSystem(grid.getNumberOfVertices())
AssemblyStiffnessMatrix(ls, grid, modulus, 0)
AssemblyGravityToVector(ls, grid, density, g, 0)
ls.applyDirichlet(0, 0)
sigma = 1e4
ls.applyNeumann(-1, -sigma)
print grid.getNumberOfVertices()
print ls.getMatrix()
print ls.getVector()
ls.solve()
print ls.getSolution()
# -----------------------------------------------------

# ------------- ANALYTICAL SOLUTION -------------------
def analyticalSolution(M, stress, L, x, gravity, rho):
	x = np.array(x)
	return x*(stress + rho*g*L)/M - rho*g*x*x/(2*M)
x_a = np.linspace(0, L, 100)
u_a = analyticalSolution(M, sigma, L, x_a, g, rho)*1000
# -----------------------------------------------------

# -------------- PLOT SOLUTION ------------------------
x_n = [v.getCoordinate() for v in grid.getVertices()]
u_n = ls.getSolution()*1000
pl.plot(u_n, x_n, 'o', label='Numeric')
pl.plot(u_a, x_a, '-', label='Analytic')
pl.grid(True)
pl.xlabel('u')
pl.ylabel('x')
pl.show()
# -----------------------------------------------------