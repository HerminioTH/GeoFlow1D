from GridLib import *
from FieldsLib import *
from FlowLib import *
from GeoLib import *
from LinearSystemLib import *
from CycleControllersLib import *
from ResultsHandlerLib import *
from UtilitiesLib import getJsonData, plotMatrix


def computeUndrainedSolution(grid, folder_settings):
	n = grid.getNumberOfVertices()

	# -------------- NUMERICAL SETTINGS -------------------
	num_set = getJsonData(folder_settings + "numerical_settings.json")
	initialTime = num_set.get("TransientCycle").get("InitialTime")
	timeStep = num_set.get("TransientCycle").get("FinalTime")
	finalTime = 10*timeStep
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
	ic = getJsonData(folder_settings + "IC.json")
	p_old = ScalarField(n, ic.get("Initial Condition").get("Value_u"))
	u_old = ScalarField(n, ic.get("Initial Condition").get("Value_p"))
	g = ic.get("Gravity")
	# -----------------------------------------------------

	# ------------- CREATE LINEAR SYSTEM ------------------
	ls = LinearSystem(2*grid.getNumberOfVertices())
	pShift = 0
	uShift = (1-pShift)
	# -----------------------------------------------------

	# ------------- BOUNDARY CONDITIONS -------------------
	bound_p = getJsonData(folder_settings + "BC_p.json")
	bound_u = getJsonData(folder_settings + "BC_u.json")
	# -----------------------------------------------------

	# --------------- FLUID FLOW MODEL --------------------
	AssemblyDarcyVelocitiesToMatrix(ls, grid, mu, k, pShift)
	AssemblyBiotAccumulationToMatrix(ls, grid, timeStep, biot, phi, cs, cf, pShift)
	AssemblyVolumetricStrainToMatrix(ls, grid, timeStep, biot, pShift)
	# ls.applyBoundaryConditionsToMatrix(grid, bound_p, pShift)
	# -----------------------------------------------------

	# -------------- GEOMECHANICAL MODEL ------------------
	AssemblyStiffnessMatrix(ls, grid, M, uShift)
	AssemblyPorePressureToMatrix(ls, grid, biot, uShift)
	ls.applyBoundaryConditionsToMatrix(grid, bound_u, uShift)
	# -----------------------------------------------------


	# -------------- TRANSIENT SOLUTION -------------------
	timeHandler.advanceTime()
	while timeHandler.isFinalTimeReached():
		ls.eraseVector()

		# --------------- FLUID FLOW MODEL --------------------
		AssemblyDarcyVelocitiesToVector(ls, grid, mu, k, rho_f, g, pShift)
		AssemblyBiotAccumulationToVector(ls, grid, timeStep, biot, phi, cs, cf, p_old, pShift)
		AssemblyVolumetricStrainToVector(ls, grid, timeStep, biot, u_old, pShift)
		# ls.applyBoundaryConditionsToVector(grid, bound_p, pShift)
		# -----------------------------------------------------

		# -------------- GEOMECHANICAL MODEL ------------------
		AssemblyGravityToVector(ls, grid, rho, g, uShift)
		ls.applyBoundaryConditionsToVector(grid, bound_u, uShift)
		# -----------------------------------------------------

		ls.solve()

		if pShift == 0:	p_new, u_new = ls.splitSolution(n)
		else:			u_new, p_new = ls.splitSolution(n)

		p_old.setField(p_new)
		u_old.setField(u_new)

		timeHandler.advanceTime()

	return p_old, u_old
	# -----------------------------------------------------

if __name__ == "__main__":
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

	folder_settings = "settings\\"
	p_old, u_old = computeUndrainedSolution(grid, folder_settings)