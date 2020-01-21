from GridLib import *
from FieldsLib import *
from LinearSystem import *
from FlowLib import *
from UtilitiesLib import *
import numpy as np
import pylab as pl

# -------------- GRID DATA ----------------------------
L_0 = 4.
L_1 = 6.
L = L_0 + L_1
nVertices = 8
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
q = -9.e1
# -----------------------------------------------------


# - BEGIN ------ NUMERICAL SOLUTION -------------------
ls = LinearSystem(grid.getNumberOfVertices())
AssemblyMassDarcyVelocities(ls, grid, mu, permeability, rho, g, pShift=0)
p_bar = 2.
Flux = 1e1

# ----------- Neumann ------------
ls.applyNeumann(-1, Flux)
for v in grid.getVertices():
	ls.applyNeumann(v.getIndex(), q*v.getVolume())
# ----- End - Neumann ------------

# ---------- Dirichlet -----------
ls.applyDirichlet(0, p_bar)
# ---- End - Dirichlet -----------

ls.solve()
# - END -------- NUMERICAL SOLUTION -------------------


# ------------- ANALYTICAL SOLUTION -------------------
def analyticalSolution(x, k, mu, q, Flux, p_bar):
	x = np.array(x)
	return (-mu*q/k/2.)*x*x + (Flux + q*L)*x*mu/k + p_bar

x_a = np.linspace(0, L, 100)
p_a = analyticalSolution(x_a, K, mu, q, Flux, p_bar)
# -----------------------------------------------------

# -------------- PLOT SOLUTION ------------------------
x_n = [v.getCoordinate() for v in grid.getVertices()]
p_n = ls.getSolution()
pl.plot(p_n, x_n, 'o', label='Numeric')
pl.plot(p_a, x_a, '-', label='Analytic')
pl.grid(True)
pl.xlabel('pressure')
pl.ylabel('x')
pl.show()
# -----------------------------------------------------