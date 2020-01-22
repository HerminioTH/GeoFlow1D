from GridLib import *
from FieldsLib import *
from FlowLib import *
from LinearSystemLib import *
from ParserLib import *
from UtilitiesLib import getJsonData

# -------------- GRID DATA ----------------------------
L = 10
nVertices = 4
nodesCoord, elemConn = createGridData(L, nVertices)
gridData = GridData()
gridData.setElementConnectivity(elemConn)
gridData.setNodeCoordinates(nodesCoord)
grid = Grid_1D(gridData)
# -----------------------------------------------------

# -------------- NUMERICAL SETTINGS -------------------
num_set = getJsonData("numerical_settings.json")
timeStep = num_set.get("TransientCycle").get("TimeStep")
finalTime = num_set.get("TransientCycle").get("FinalTime")
# -----------------------------------------------------

# -------------- PROPERTIES ---------------------------
fluid = getJsonData("fluid.json")
rho = fluid.get("WATER").get("Density").get("value")
mu = fluid.get("WATER").get("Viscosity").get("value")
cf = fluid.get("WATER").get("Compressibility").get("value")
g = -9.81

solid = getJsonData("solid.json")
k = ScalarField(grid.getNumberOfRegions())
phi = ScalarField(grid.getNumberOfRegions())
cs = ScalarField(grid.getNumberOfRegions())
for region in grid.getRegions():
	k.setValue(region, solid.get("ROCK_1").get("Permeability").get("value"))
	phi.setValue(region, solid.get("ROCK_1").get("Porosity").get("value"))
	cs.setValue(region, solid.get("ROCK_1").get("Compressibility").get("value"))
# -----------------------------------------------------


# -------------- NUMERICAL SOLUTION -------------------
ls = LinearSystem(grid.getNumberOfVertices())
p_init = 2.0
p_old = ScalarField(grid.getNumberOfVertices(), p_init)
AssemblyDarcyVelocitiesToMatrix(ls, grid, mu, k)
AssemblyDarcyVelocitiesToVector(ls, grid, mu, k, rho, g)
AssemblyFluidFlowAccumulationToMatrix(ls, grid, timeStep, phi, cs, cf)
AssemblyFluidFlowAccumulationToVector(ls, grid, timeStep, phi, cs, cf, p_old)

flux = -2.5
ls.applyNeumann(-1, flux)

p_bar = 1.0
ls.applyDirichlet(0, p_bar)
ls.solve()
print ls.getMatrix()
print ls.getVector()
print ls.getSolution()
# -----------------------------------------------------