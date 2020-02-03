from GridLib import *
from FieldsLib import *
from FlowLib import *
from GeoLib import *
from LinearSystemLib import *
from CycleControllersLib import *
from ResultsHandlerLib import *
from undrainedSolution import *
from UtilitiesLib import *
import numpy as np



# ------------------ GRID DATA ------------------------
L = 10
nVertices = 15
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
maxIte = num_set.get("IterativeCycle").get("MaximumNumberOfIterations")
maxTol = num_set.get("IterativeCycle").get("Tolerance")
timeHandler = TimeHandler(timeStep, finalTime, initialTime)
iterativeController = IterativeCycleController(maxIte, maxTol)
method_split = num_set.get("SplitMethod").get("Name")
folder_results += method_split + "_DELTA\\"
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
K = ScalarField(grid.getNumberOfRegions())
biot = ScalarField(grid.getNumberOfRegions())
rho_s = ScalarField(grid.getNumberOfRegions())
for region in grid.getRegions():
	k.setValue(region, solid.get("ROCK_1").get("Permeability").get("value"))
	phi.setValue(region, solid.get("ROCK_1").get("Porosity").get("value"))
	cs.setValue(region, solid.get("ROCK_1").get("Compressibility").get("value"))
	rho_s.setValue(region, solid.get("ROCK_1").get("Density").get("value"))
	G = solid.get("ROCK_1").get("ShearModulus").get("value")
	nu = solid.get("ROCK_1").get("PoissonsRatio").get("value")
	bulk = 2*G*(1 + nu)/(3*(1 - 2*nu))
	M.setValue(region, bulk + 4*G/3.)
	K.setValue(region, bulk)
	CS = solid.get("ROCK_1").get("Compressibility").get("value")
	biot.setValue(region, 1 - CS*bulk)

rho = ScalarField(grid.getNumberOfRegions())
for r in grid.getRegions():
	rho.setValue(r, phi.getValue(r)*rho_f + (1 - phi.getValue(r))*rho_s.getValue(r))
delta = ScalarField(grid.getNumberOfVertices(), num_set.get("SplitMethod").get("Relaxation"))
# delta = ScalarField(grid.getNumberOfVertices(), 1. + 4*G/bulk/3.)
# delta = ScalarField(grid.getNumberOfVertices(), 1.0)
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
res_error = SaveResults(grid, "error.txt", folder_results, 'Error', '-')
res_u.copySettings(folder_settings, folder_results)
# -----------------------------------------------------

# -------------- TRANSIENT SOLUTION -------------------
res_p.saveField(timeHandler.getCurrentTime(), p_old.getField())
res_u.saveField(timeHandler.getCurrentTime(), u_old.getField())
timeHandler.advanceTime()

rates = []
isFirst = True
while timeHandler.isFinalTimeReached():
##	timeHandler.printCurrentTime()
	iterativeController.reset()
	error_list = []
	while iterativeController.keepCycling():

		# --------------- FLUID FLOW MODEL --------------------
		ls_mass.resetLinearSystem()

		AssemblyBiotAccumulationToMatrix(ls_mass, grid, timeStep, biot, phi, cs, cf, pShift)
		AssemblyBiotAccumulationToVector(ls_mass, grid, timeStep, biot, phi, cs, cf, p_old, pShift)

		if method_split == "FIXED_STRESS":
			AssemblyFixedStressAccumulationToMatrix(ls_mass, grid, timeStep, biot, delta, K, pShift)
			AssemblyFixedStressAccumulationToVector(ls_mass, grid, timeStep, biot, delta, K, p_new, pShift)

		AssemblyDarcyVelocitiesToMatrix(ls_mass, grid, mu, k, pShift)
		AssemblyDarcyVelocitiesToVector(ls_mass, grid, mu, k, rho_f, g, pShift)

		AssemblyVolumetricStrainToVector(ls_mass, grid, timeStep, biot, u_old, pShift)
		AssemblyVolumetricStrainToVector(ls_mass, grid, -timeStep, biot, u_new, pShift)

		ls_mass.applyBoundaryConditionsToMatrix(grid, bound_p, pShift)
		ls_mass.applyBoundaryConditionsToVector(grid, bound_p, pShift)

		ls_mass.solve()
		error = p_new.getField() - ls_mass.getSolution()
		L2_mass = computeNormL2(error, grid)
		p_new.setField(ls_mass.getSolution())
		# -----------------------------------------------------


		# -------------- GEOMECHANICAL MODEL ------------------
		ls_geom.resetLinearSystem()
		AssemblyStiffnessMatrix(ls_geom, grid, M, uShift)
		AssemblyGravityToVector(ls_geom, grid, rho, g, uShift)
		AssemblyPorePressureToVector(ls_geom, grid, biot, p_new, uShift)
		ls_geom.applyBoundaryConditionsToMatrix(grid, bound_u, uShift)
		ls_geom.applyBoundaryConditionsToVector(grid, bound_u, uShift)

		ls_geom.solve()
		error = u_new.getField() - ls_geom.getSolution()
		L2_geom = computeNormL2(error, grid)
		u_new.setField(ls_geom.getSolution())
		# -----------------------------------------------------

		error_list.append(L2_mass)
		iterativeController.execute(max(L2_mass, L2_geom))

	r = (np.log10(error_list[1]) - np.log10(error_list[-1]))/len(error_list)
	rates.append(r)
	# print error_list
	try:		m5 = computeMedia(rates, 5)
	except:		m5 = r

	try:		m10 = computeMedia(rates, 10)
	except:		m10 = r
	print timeHandler.getCurrentTime(), len(error_list), r, m5, m10


	if not isFirst:
		d = float(raw_input("Delta: "))
		delta = ScalarField(grid.getNumberOfVertices(), d)
	isFirst = False

	res_error.saveField(timeHandler.getCurrentTime(), np.array(error_list))

##	print "Ite: %i"%iterativeController.iteNumber
##	print "L2_mass: %e"%L2_mass
##	print "L2_geom: %e"%L2_geom


	p_old.setField(p_new.getField())
	u_old.setField(u_new.getField())

	res_p.saveField(timeHandler.getCurrentTime(), p_old.getField())
	res_u.saveField(timeHandler.getCurrentTime(), u_old.getField())

	timeHandler.advanceTime()

res_p.close()
res_u.close()
# -----------------------------------------------------

# print ls.getMatrix()
# print ls.getVector()

print "P_eq num = ",p_new.getField()[0]

top_stress = getJsonData(folder_settings + "BC_u.json").get("TOP").get("Value")
for region in grid.getRegions():
	Q = 1/( cs.getValue(region)*(biot.getValue(region) - phi.getValue(region)) + cf*phi.getValue(region) )
	p_eq = biot.getValue(region)*Q*top_stress/(M.getValue(region) + Q*biot.getValue(region)**2)
print "Equilibrium pressure is: %f Pa"%p_eq

print "Delta = ", 1. + 4*G/bulk/3.
