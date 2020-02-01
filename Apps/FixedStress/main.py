from GridLib import *
from FieldsLib import *
from FlowLib import *
from GeoLib import *
from LinearSystemLib import *
from CycleControllersLib import *
from ResultsHandlerLib import *
from undrainedSolution import *
from UtilitiesLib import *

# ------------------ GRID DATA ------------------------
L = 10
nVertices = 7
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
folder_settings = "settings\\"
folder_results = "results\\"
# -----------------------------------------------------

# -------------- NUMERICAL SETTINGS -------------------
num_set = getJsonData(folder_settings + "numerical_settings.json")
initialTime = num_set.get("TransientCycle").get("InitialTime")
timeStep = num_set.get("TransientCycle").get("TimeStep")
finalTime = num_set.get("TransientCycle").get("FinalTime")
timeHandler = TimeHandler(timeStep, finalTime, initialTime)
iterativeController = IterativeCycleController(500, 1e-1)
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
p_old = ScalarField(grid.getNumberOfVertices())
u_old = ScalarField(grid.getNumberOfVertices())
p_new = ScalarField(grid.getNumberOfVertices())
u_new = ScalarField(grid.getNumberOfVertices())
ic = getJsonData(folder_settings + "IC.json")
g = ic.get("Gravity")
# -----------------------------------------------------

# ------------- CREATE LINEAR SYSTEM ------------------
ls_mass = LinearSystem(grid.getNumberOfVertices())
ls_geom = LinearSystem(grid.getNumberOfVertices())
pShift = uShift = 0
n = grid.getNumberOfVertices()
# -----------------------------------------------------

# ------------- BOUNDARY CONDITIONS -------------------
bound_p = getJsonData(folder_settings + "BC_p.json")
bound_u = getJsonData(folder_settings + "BC_u.json")
# -----------------------------------------------------




# --------------- RESULTS HANDLER ---------------------
res_p = SaveResults(grid, "p.txt", folder_results, 'Pressure', 'Pa')
res_u = SaveResults(grid, "u.txt", folder_results, 'Displacement', 'm')
# res_u.copySettings(folder_settings, folder_results)
# -----------------------------------------------------

# -------------- TRANSIENT SOLUTION -------------------
res_p.saveField(timeHandler.getCurrentTime(), p_old.getField())
res_u.saveField(timeHandler.getCurrentTime(), u_old.getField())
timeHandler.advanceTime()

while timeHandler.isFinalTimeReached():
	timeHandler.printCurrentTime()

	ls_mass.resetLinearSystem()
	ls_geom.resetLinearSystem()
	iterativeController.reset()
	while iterativeController.keepCycling():
		print u_old.getField()
		print u_new.getField()
		# --------------- FLUID FLOW MODEL --------------------
		AssemblyDarcyVelocitiesToMatrix(ls_mass, grid, mu, k, pShift)
		AssemblyBiotAccumulationToMatrix(ls_mass, grid, timeStep, biot, phi, cs, cf, pShift)
		AssemblyDarcyVelocitiesToVector(ls_mass, grid, mu, k, rho_f, g, pShift)
		AssemblyBiotAccumulationToVector(ls_mass, grid, timeStep, biot, phi, cs, cf, p_old, pShift)
		AssemblyVolumetricStrainToVector(ls_mass, grid, timeStep, biot, u_old, pShift)
		AssemblyVolumetricStrainToVector(ls_mass, grid, -timeStep, biot, u_new, pShift)
		ls_mass.applyBoundaryConditionsToMatrix(grid, bound_p, pShift)
		ls_mass.applyBoundaryConditionsToVector(grid, bound_p, pShift)


		ls_mass.solve()
		error = p_new.getField() - ls_mass.getSolution()
		L2_mass = computeNormL2(error, grid)
		p_new.setField(ls_mass.getSolution())

		# print L2_mass
		print ls_mass.getVector()
		# -----------------------------------------------------

		# -------------- GEOMECHANICAL MODEL ------------------
		AssemblyStiffnessMatrix(ls_geom, grid, M, uShift)
		AssemblyGravityToVector(ls_geom, grid, rho, g, uShift)
		ls_geom.applyBoundaryConditionsToMatrix(grid, bound_u, uShift)
		ls_geom.applyBoundaryConditionsToVector(grid, bound_u, uShift)

		ls_geom.solve()
		u_new.setField(ls_geom.getSolution())
		# -----------------------------------------------------

		iterativeController.execute(L2_mass)

	print '-------------------'
	print p_old.getField()
	print p_new.getField()
	print '-------------------'
	print u_old.getField()
	print u_new.getField(),'\n'

	p_old.setField(p_new.getField())
	u_old.setField(u_new.getField())



	res_p.saveField(timeHandler.getCurrentTime(), p_old.getField())
	res_u.saveField(timeHandler.getCurrentTime(), u_old.getField())

	# p_old.setField(p_new)
	# u_old.setField(u_new)

	timeHandler.advanceTime()

res_p.close()
res_u.close()
# -----------------------------------------------------

# print ls.getMatrix()
# print ls.getVector()