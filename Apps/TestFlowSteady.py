from GridLib import *
from FieldsLib import *
from LinearSystem import *
from FlowLib import *
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
K = 1.2e-12
permeability = ScalarField(grid.getNumberOfRegions())
permeability.setValue(grid.getRegions()[0], K)
rho = 1000.
mu =1e-3
g = -9.81
# -----------------------------------------------------

# -------------- NUMERICAL SOLUTION -------------------
ls = LinearSystem(grid.getNumberOfVertices())
AssemblyMassDarcyVelocities(ls, grid, mu, permeability, density=rho, gravity=g, pShift=0)
# ls.applyDirichlet(0, 0)
# ls.applyDirichlet(-1, 10)
# ls.applyDirichlet(0, 0)
ls.applyNeumann(0, -9e-12)
ls.applyNeumann(-1, -9e-12)
ls.applyDirichlet(-8, 0)
ls.solve()
# -----------------------------------------------------

# # ------------- ANALYTICAL SOLUTION -------------------
# def analyticalSolution(M, stress, L, x, gravity, rho):
# 	x = np.array(x)
# 	return x*(stress + rho*g*L)/M - rho*g*x*x/(2*M)
# x_a = np.linspace(0, L, 100)
# u_a = analyticalSolution(M, sigma, L, x_a, g, rho)*1000
# # -----------------------------------------------------

# -------------- PLOT SOLUTION ------------------------
x_n = [v.getCoordinate() for v in grid.getVertices()]
u_n = ls.getSolution()*1000
pl.plot(u_n, x_n, 'o', label='Numeric')
# pl.plot(u_a, x_a, '-', label='Analytic')
pl.grid(True)
pl.xlabel('u')
pl.ylabel('x')
pl.show()
# -----------------------------------------------------