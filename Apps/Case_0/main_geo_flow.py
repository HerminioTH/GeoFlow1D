from GridLib import *
from FieldsLib import *
from FlowLib import *
from GeoLib import *
from LinearSystemLib import *
from CycleControllersLib import *
from ResultsHandlerLib import *
from UtilitiesLib import getJsonData, plotMatrix

# ------------------ GRID DATA ------------------------
L = 10
nVertices = 60
nodesCoord, elemConn = createGridData(L, nVertices)
gridData = GridData()
gridData.setElementConnectivity(elemConn)
gridData.setNodeCoordinates(nodesCoord)
gridData.addBoundary("TOP", nVertices-2, nVertices-1)
gridData.addBoundary("BOTTOM", 0, 0)
# -----------------------------------------------------

# --------------------- GRID --------------------------
grid = Grid_1D(gridData)
# -----------------------------------------------------

# ---------------- FOLDER SETTINGS --------------------
folder_settings = "settings/"
folder_results = "results/geo_flow/"
# -----------------------------------------------------

# -------------- NUMERICAL SETTINGS -------------------
num_set = getJsonData(folder_settings + "numerical_settings.json")
initialTime = num_set.get("TransientCycle").get("InitialTime")
timeStep = num_set.get("TransientCycle").get("TimeStep")
finalTime = num_set.get("TransientCycle").get("FinalTime")
timeHandler = TimeHandler( timeStep, finalTime, initialTime )
# -----------------------------------------------------

# -------------- PROPERTIES ---------------------------
fluid = getJsonData(folder_settings + "properties_fluid.json")
rho_f = fluid.get("WATER").get("Density").get("value")
mu = fluid.get("WATER").get("Viscosity").get("value")
cf = fluid.get("WATER").get("Compressibility").get("value")

solid = getJsonData(folder_settings + "properties_solid.json")
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

rho = ScalarField(grid.getNumberOfRegions())
for r in grid.getRegions():
	rho.setValue(r, phi.getValue(r)*rho_f + (1 - phi.getValue(r))*rho_s.getValue(r))
# -----------------------------------------------------

# ---------------- INITIAL FIELDS ---------------------
p_old = ScalarField(grid.getNumberOfVertices(), 0.0)
u_old = ScalarField(grid.getNumberOfVertices(), 0.0)
ic = getJsonData(folder_settings + "IC.json")
g = ic.get("Gravity")
# -----------------------------------------------------

# ------------- CREATE LINEAR SYSTEM ------------------
ls = LinearSystem(2*grid.getNumberOfVertices())
pShift = 1
uShift = (1-pShift)
n = grid.getNumberOfVertices()
# -----------------------------------------------------

# ------------- BOUNDARY CONDITIONS -------------------
bound_p = getJsonData(folder_settings + "BC_p.json")
bound_u = getJsonData(folder_settings + "BC_u.json")
# -----------------------------------------------------

# --------------- FLUID FLOW MODEL --------------------
AssemblyDarcyVelocitiesToMatrix(ls, grid, mu, k, pShift)
AssemblyBiotAccumulationToMatrix(ls, grid, timeStep, biot, phi, cs, cf, pShift)
AssemblyVolumetricStrainToMatrix(ls, grid, timeStep, biot, pShift)
ls.applyBoundaryConditionsToMatrix(grid, bound_p, pShift)
# -----------------------------------------------------

# -------------- GEOMECHANICAL MODEL ------------------
AssemblyStiffnessMatrix(ls, grid, M, uShift)
AssemblyPorePressureToMatrix(ls, grid, biot, uShift)
ls.applyBoundaryConditionsToMatrix(grid, bound_u, uShift)
# -----------------------------------------------------

# --------------- RESULTS HANDLER ---------------------
res_p = SaveResults(grid, "p.txt", folder_results, 'Pressure', 'Pa')
res_u = SaveResults(grid, "u.txt", folder_results, 'Displacement', 'm')
res_u.copySettings(folder_settings, folder_results)
# -----------------------------------------------------

# -------------- TRANSIENT SOLUTION -------------------
res_p.saveField(timeHandler.getCurrentTime(), p_old.getField())
res_u.saveField(timeHandler.getCurrentTime(), u_old.getField())
timeHandler.advanceTime()
while timeHandler.isFinalTimeReached():
	print(timeHandler.getCurrentTime())
	ls.eraseVector()

	# --------------- FLUID FLOW MODEL --------------------
	AssemblyDarcyVelocitiesToVector(ls, grid, mu, k, rho_f, g, pShift)
	AssemblyBiotAccumulationToVector(ls, grid, timeStep, biot, phi, cs, cf, p_old, pShift)
	AssemblyVolumetricStrainToVector(ls, grid, timeStep, biot, u_old, pShift)
	ls.applyBoundaryConditionsToVector(grid, bound_p, pShift)
	# -----------------------------------------------------

	# -------------- GEOMECHANICAL MODEL ------------------
	AssemblyGravityToVector(ls, grid, rho, g, uShift)
	ls.applyBoundaryConditionsToVector(grid, bound_u, uShift)
	# -----------------------------------------------------

	ls.solve()

	if pShift == 0:	p_new, u_new = ls.splitSolution(n)
	else:			u_new, p_new = ls.splitSolution(n)

	res_p.saveField(timeHandler.getCurrentTime(), p_new)
	res_u.saveField(timeHandler.getCurrentTime(), u_new)

	p_old.setField(p_new)
	u_old.setField(u_new)

	timeHandler.advanceTime()

res_p.close()
res_u.close()
# -----------------------------------------------------

# print ls.getMatrix()
# print ls.getVector()
# plotMatrix(ls.getMatrix())
