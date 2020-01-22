from GridLib import *
from FieldsLib import *
from FlowLib import *
from GeoLib import *
from LinearSystemLib import *
from ParserLib import *
from TimeHandlerLib import *
from ResultsHandlerLib import *
from UtilitiesLib import getJsonData, plotMatrix

# -------------- GRID DATA ----------------------------
L = 10
nVertices = 40
nodesCoord, elemConn = createGridData(L, nVertices)
gridData = GridData()
gridData.setElementConnectivity(elemConn)
gridData.setNodeCoordinates(nodesCoord)
grid = Grid_1D(gridData)
# -----------------------------------------------------


# -------------- NUMERICAL SETTINGS -------------------
num_set = getJsonData("settings\\numerical_settings.json")
initialTime = num_set.get("TransientCycle").get("InitialTime")
timeStep = num_set.get("TransientCycle").get("TimeStep")
finalTime = num_set.get("TransientCycle").get("FinalTime")
timeHandler = TimeHandler( timeStep, finalTime, initialTime )
# -----------------------------------------------------

# -------------- PROPERTIES ---------------------------
fluid = getJsonData("properties\\fluid.json")
rho_f = fluid.get("WATER").get("Density").get("value")
mu = fluid.get("WATER").get("Viscosity").get("value")
cf = fluid.get("WATER").get("Compressibility").get("value")
g = 0.0

solid = getJsonData("properties\\solid.json")
k = ScalarField(grid.getNumberOfRegions())
phi = ScalarField(grid.getNumberOfRegions())
cs = ScalarField(grid.getNumberOfRegions())
M = ScalarField(grid.getNumberOfRegions())
biot = ScalarField(grid.getNumberOfRegions())
rho_s = ScalarField(grid.getNumberOfRegions())
for region in grid.getRegions():
	k.setValue(region, solid.get("ROCK_1").get("Permeability").get("value"))
	phi.setValue(region, solid.get("ROCK_1").get("Porosity").get("value"))
	cs.setValue(region, solid.get("ROCK_1").get("Compressibility").get("value"))
	rho_s.setValue(region, solid.get("ROCK_1").get("Density").get("value"))
	G = solid.get("ROCK_1").get("ShearModulus").get("value")
	nu = solid.get("ROCK_1").get("PoissonsRatio").get("value")
	K = 2*G*(1 + nu)/(3*(1 - 2*nu))
	M.setValue(region, K + 4*G/3.)
	CS = solid.get("ROCK_1").get("Compressibility").get("value")
	biot.setValue(region, 1 - CS*K)
# -----------------------------------------------------

# --------------- RESULTS HANDLER ---------------------
folderName = "results\\"
res_p = SaveResults(grid, "p.txt", folderName, 'Pressure', 'Pa')
res_u = SaveResults(grid, "u.txt", folderName, 'Displacement', 'm')
# -----------------------------------------------------

# ---------------- INITIAL FIELDS ---------------------
p_old = ScalarField(grid.getNumberOfVertices(), 50e3)
u_old = ScalarField(grid.getNumberOfVertices(), 0.0)
# -----------------------------------------------------

# ------------- CREATE LINEAR SYSTEM ------------------
ls = LinearSystem(2*grid.getNumberOfVertices())
pShift = 1
uShift = 0
n = grid.getNumberOfVertices()
# -----------------------------------------------------

# --------------- FLUID FLOW MODEL --------------------
AssemblyDarcyVelocitiesToMatrix(ls, grid, mu, k, pShift)
AssemblyDarcyVelocitiesToVector(ls, grid, mu, k, rho_f, g, pShift)
AssemblyFluidFlowAccumulationToMatrix(ls, grid, timeStep, phi, cs, cf, pShift)
AssemblyFluidFlowAccumulationToVector(ls, grid, timeStep, phi, cs, cf, p_old, pShift)
p_top = 0.0
ls.applyDirichlet(2*n-1, p_top)
# -----------------------------------------------------

# -------------- GEOMECHANICAL MODEL ------------------
AssemblyStiffnessMatrix(ls, grid, M, uShift)
AssemblyGravityToVector(ls, grid, rho_s, g, uShift)
AssemblyPorePressureToMatrix(ls, grid, biot, uShift)
u_bottom = 0.0
top_stress = 1.0e3
ls.applyNeumann(n-1, top_stress)
ls.applyDirichlet(0, u_bottom)
# -----------------------------------------------------



# -------------- TRANSIENT SOLUTION -------------------
timeHandler.advanceTime()
while timeHandler.isFinalTimeReached():
	ls.solve()
	x_u, x_p = ls.splitSolution(n)
	p_old.setField(x_p)
	u_old.setField(x_u)
	res_p.saveField(timeHandler.getCurrentTime(), p_old.getField())
	res_u.saveField(timeHandler.getCurrentTime(), u_old.getField())

	ls.eraseVector()

	# --------------- FLUID FLOW MODEL --------------------
	AssemblyDarcyVelocitiesToVector(ls, grid, mu, k, rho_f, g, pShift)
	AssemblyFluidFlowAccumulationToVector(ls, grid, timeStep, phi, cs, cf, p_old, pShift)
	ls.applyDirichletToVector(2*n-1, p_top)
	# -----------------------------------------------------

	# -------------- GEOMECHANICAL MODEL ------------------
	AssemblyGravityToVector(ls, grid, rho_s, g, uShift)
	ls.applyNeumann(n-1, top_stress)
	ls.applyDirichletToVector(0, u_bottom)
	# -----------------------------------------------------

	timeHandler.advanceTime()
res_p.close()
res_u.close()
# -----------------------------------------------------

# print ls.getMatrix()
# print ls.getVector()
plotMatrix(ls.getMatrix())