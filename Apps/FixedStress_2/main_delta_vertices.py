from GridLib import *
from FieldsLib import *
from FlowLib import *
from GeoLib import *
from LinearSystemLib import *
from CycleControllersLib import *
from ResultsHandlerLib import *
from undrainedSolution import *
from UtilitiesLib import *
from Optimization import *
import numpy as np



def computeDimensionalTime(fluid, solid, L):
	mu = fluid.get("WATER").get("Viscosity").get("value")
	CF = fluid.get("WATER").get("Compressibility").get("value")
	for region in grid.getRegions():
		G = solid.get("ROCK_1").get("ShearModulus").get("value")
		nu = solid.get("ROCK_1").get("PoissonsRatio").get("value")
		bulk = 2*G*(1 + nu)/(3*(1 - 2*nu))
		CS = solid.get("ROCK_1").get("Compressibility").get("value")
		phi_value = solid.get("ROCK_1").get("Porosity").get("value")
		alpha = 1 - CS*bulk
		mv = 1.0/(bulk + (4.0/3.0)*G)
		k = solid.get("ROCK_1").get("Permeability").get("value")
		Q = 1.0/(phi_value*CF + CS*(alpha - phi_value))
		c_v = (Q/((alpha**2)*Q*mv + 1.0))*(k/mu)
		time = (L**2)/c_v
	return time



# ------------------ GRID DATA ------------------------
L = 10
nVertices = 6
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
folder_results = "results\\TESTE\\"
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
	CS = solid.get("ROCK_1").get("Compressibility").get("value")
	bulk = 2*G*(1 + nu)/(3*(1 - 2*nu))
	pWave = bulk + 4*G/3.
	alpha = 1 - CS*bulk
	M.setValue(region, pWave)
	K.setValue(region, bulk)
	biot.setValue(region, alpha)
rho = ScalarField(grid.getNumberOfRegions())
for r in grid.getRegions():
	rho.setValue(r, phi.getValue(r)*rho_f + (1 - phi.getValue(r))*rho_s.getValue(r))
# -----------------------------------------------------

# -------------- NUMERICAL SETTINGS -------------------
num_set = getJsonData(folder_settings + "numerical_settings.json")
initialTime = 0.0
finalTime = computeDimensionalTime(fluid, solid, L)
timeStep = finalTime/500.
maxIte = num_set.get("IterativeCycle").get("MaximumNumberOfIterations")
maxTol = num_set.get("IterativeCycle").get("Tolerance")
timeHandler = TimeHandler(timeStep, finalTime, initialTime)
iterativeController = IterativeCycleController(maxIte, maxTol)
# -----------------------------------------------------

# -------------- SPLITTING SETTINGS -------------------
step = num_set.get("SplitMethod").get("Step")
reductionFactor = num_set.get("SplitMethod").get("Factor")
mediaMovel1 = num_set.get("SplitMethod").get("Media1")
mediaMovel2 = num_set.get("SplitMethod").get("Media2")
d = num_set.get("SplitMethod").get("Relaxation")
delta = BackForwardStep(d, step, reductionFactor)
deltaField = ScalarField(grid.getNumberOfVertices(), d)
method_split = num_set.get("SplitMethod").get("Name")
manual = num_set.get("SplitMethod").get("Manual")
if method_split == "FIXED_STRESS_D":
	folder_results += method_split
	folder_results += "_" + str(step)
	folder_results += "_" + str(reductionFactor)
	folder_results += "_" + str(mediaMovel1)
	folder_results += "_" + str(mediaMovel2)
	folder_results += "\\"
else:
	folder_results += method_split + "\\"
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
res_delta = SaveResults(grid, "delta.txt", folder_results, 'delta', '-')
res_u.copySettings(folder_settings, folder_results)
# -----------------------------------------------------

# -------------- TRANSIENT SOLUTION -------------------
res_p.saveField(timeHandler.getCurrentTime(), p_old.getField())
res_u.saveField(timeHandler.getCurrentTime(), u_old.getField())
res_delta.saveField(timeHandler.getCurrentTime(), deltaField.getField())
timeHandler.advanceTime()

rates = []
rate_vertices_old = np.zeros(grid.getNumberOfVertices())
counter = 0
while timeHandler.isFinalTimeReached():
##	timeHandler.printCurrentTime()
	iterativeController.reset()
	error_list = []
	while iterativeController.keepCycling():

		# --------------- FLUID FLOW MODEL --------------------
		ls_mass.resetLinearSystem()
		AssemblyBiotAccumulationToMatrix(ls_mass, grid, timeStep, biot, phi, cs, cf, pShift)
		AssemblyBiotAccumulationToVector(ls_mass, grid, timeStep, biot, phi, cs, cf, p_old, pShift)
		if method_split == "FIXED_STRESS_K" or method_split == "FIXED_STRESS_D":
			AssemblyFixedStressAccumulationToMatrix(ls_mass, grid, timeStep, biot, deltaField, K, pShift)
			AssemblyFixedStressAccumulationToVector(ls_mass, grid, timeStep, biot, deltaField, K, p_new, pShift)
		elif method_split == "FIXED_STRESS_M":
			AssemblyFixedStressAccumulationToMatrix(ls_mass, grid, timeStep, biot, deltaField, M, pShift)
			AssemblyFixedStressAccumulationToVector(ls_mass, grid, timeStep, biot, deltaField, M, p_new, pShift)
		else:
			pass
		AssemblyDarcyVelocitiesToMatrix(ls_mass, grid, mu, k, pShift)
		AssemblyDarcyVelocitiesToVector(ls_mass, grid, mu, k, rho_f, g, pShift)
		AssemblyVolumetricStrainToVector(ls_mass, grid, timeStep, biot, u_old, pShift)
		AssemblyVolumetricStrainToVector(ls_mass, grid, -timeStep, biot, u_new, pShift)
		ls_mass.applyBoundaryConditionsToMatrix(grid, bound_p, pShift)
		ls_mass.applyBoundaryConditionsToVector(grid, bound_p, pShift)
		ls_mass.solve()
		# -----------------------------------------------------

		# ----------------- RELATIVE ERROR --------------------
		numerator = computeNormL2(p_new.getField() - ls_mass.getSolution(), grid)
		if len(error_list) == 0:
			error_vertices_old = p_old.getField() - ls_mass.getSolution()
			denominator = computeNormL2(error_vertices_old, grid)
			if denominator == 0.: denominator = 1
		L2_mass = numerator/denominator
		# -----------------------------------------------------

		# ---------------- RELATIVE ERROR - 2 -----------------
		if iterativeController.iteNumber < 1:
			den_vertices = abs(p_old.getField() - ls_mass.getSolution())
			for i in range(den_vertices.size):
				if den_vertices[i] == 0.0:
					den_vertices[i] = 1.0
		error_vertices = abs(p_new.getField() - ls_mass.getSolution())/den_vertices
		for i in range(den_vertices.size):
			if error_vertices[i] == 0.0:
				error_vertices[i] = 1.0
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
		u_new.setField(ls_geom.getSolution())
		# -----------------------------------------------------

		error_list.append(L2_mass)
		iterativeController.execute(L2_mass)

	rate_vertices_new = np.log10(error_vertices)/iterativeController.iteNumber
	rate_rate = -(rate_vertices_new - rate_vertices_old)
	rate_vertices_old = rate_vertices_new

	rate = computeRate(error_list)
	rates.append(rate)

	print len(error_list), d, rate, rate_rate

	if method_split == "FIXED_STRESS_D" and manual == False:
		d = delta.computeDelta(rate)
		deltaField = ScalarField(grid.getNumberOfVertices(), d)





	# print timeHandler.getCurrentTime(), len(error_list), d, rate, denominator, L2_mass


	res_error.saveField(timeHandler.getCurrentTime(), np.array(error_list))
	res_delta.saveField(timeHandler.getCurrentTime(), deltaField.getField())

	p_old.setField(p_new.getField())
	u_old.setField(u_new.getField())

	res_p.saveField(timeHandler.getCurrentTime(), p_old.getField())
	res_u.saveField(timeHandler.getCurrentTime(), u_old.getField())

	timeHandler.advanceTime()
	counter += 1

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
print folder_results
