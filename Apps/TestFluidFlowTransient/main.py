from GridLib import *
from FieldsLib import *
from FlowLib import *
from LinearSystemLib import *
from CycleControllersLib import *
from ResultsHandlerLib import *
from UtilitiesLib import getJsonData

# -------------- GRID DATA ----------------------------
L = 10
nVertices = 20
nodesCoord, elemConn = createGridData(L, nVertices)
gridData = GridData()
gridData.setElementConnectivity(elemConn)
gridData.setNodeCoordinates(nodesCoord)
grid = Grid_1D(gridData)
# -----------------------------------------------------

# -------------- NUMERICAL SETTINGS -------------------
num_set = getJsonData("numerical_settings.json")
initialTime = num_set.get("TransientCycle").get("InitialTime")
timeStep = num_set.get("TransientCycle").get("TimeStep")
finalTime = num_set.get("TransientCycle").get("FinalTime")
timeHandler = TimeHandler( timeStep, finalTime, initialTime )
# -----------------------------------------------------

# -------------- PROPERTIES ---------------------------
fluid = getJsonData("fluid.json")
rho = fluid.get("WATER").get("Density").get("value")
mu = fluid.get("WATER").get("Viscosity").get("value")
cf = fluid.get("WATER").get("Compressibility").get("value")
g = 0.0

solid = getJsonData("solid.json")
k = ScalarField(grid.getNumberOfRegions())
phi = ScalarField(grid.getNumberOfRegions())
cs = ScalarField(grid.getNumberOfRegions())
for region in grid.getRegions():
	k.setValue(region, solid.get("ROCK_1").get("Permeability").get("value"))
	phi.setValue(region, solid.get("ROCK_1").get("Porosity").get("value"))
	cs.setValue(region, solid.get("ROCK_1").get("Compressibility").get("value"))
# -----------------------------------------------------

# --------------- RESULTS HANDLER ---------------------
fileName = "p.txt"
folderName = "results\\"
results = SaveResults(grid, fileName, folderName, 'Pressure', 'Pa')
# -----------------------------------------------------


# ------------ ASSEMBLY LINEAR SYSTEM -----------------
ls = LinearSystem(grid.getNumberOfVertices())
p_init = 10.0
p_old = ScalarField(grid.getNumberOfVertices(), p_init)
AssemblyDarcyVelocitiesToMatrix(ls, grid, mu, k)
AssemblyDarcyVelocitiesToVector(ls, grid, mu, k, rho, g)
AssemblyFluidFlowAccumulationToMatrix(ls, grid, timeStep, phi, cs, cf)
AssemblyFluidFlowAccumulationToVector(ls, grid, timeStep, phi, cs, cf, p_old)

p_top = 1.0
p_bottom = 5.0
ls.applyDirichlet(0, p_bottom)
# ls.applyDirichlet(-1, p_top)
ls.applyNeumann(-1, 1e-10)
# -----------------------------------------------------


# -------------- TRANSIENT SOLUTION -------------------
timeHandler.advanceTime()
while timeHandler.isFinalTimeReached():
	ls.solve()
	p_old.setField(ls.getSolution())
	results.saveField(timeHandler.getCurrentTime(), p_old.getField())

	ls.eraseVector()
	AssemblyDarcyVelocitiesToVector(ls, grid, mu, k, rho, g)
	AssemblyFluidFlowAccumulationToVector(ls, grid, timeStep, phi, cs, cf, p_old)
	ls.applyDirichletToVector(0, p_bottom)
	# ls.applyDirichletToVector(-1, p_top)
	ls.applyNeumann(-1, 1e-10)

	timeHandler.advanceTime()
results.close()
# -----------------------------------------------------