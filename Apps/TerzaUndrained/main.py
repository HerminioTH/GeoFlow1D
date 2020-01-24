from GridLib import *
from FieldsLib import *
from FlowLib import *
from GeoLib import *
from LinearSystemLib import *
from TimeHandlerLib import *
from ResultsHandlerLib import *
from UtilitiesLib import getJsonData, plotMatrix

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
num_set = getJsonData("settings\\numerical_settings.json")
initialTime = num_set.get("TransientCycle").get("InitialTime")
timeStep = num_set.get("TransientCycle").get("TimeStep")
finalTime = num_set.get("TransientCycle").get("FinalTime")
timeHandler = TimeHandler( timeStep, finalTime, initialTime )
# -----------------------------------------------------

# -------------- PROPERTIES ---------------------------
fluid = getJsonData("settings\\fluid.json")
rho_f = fluid.get("WATER").get("Density").get("value")
mu = fluid.get("WATER").get("Viscosity").get("value")
cf = fluid.get("WATER").get("Compressibility").get("value")
g = 0.0

solid = getJsonData("settings\\solid.json")
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
destinationFolder = "results\\"
sourceFolder = "settings\\"
res_p = SaveResults(grid, "p.txt", destinationFolder, 'Pressure', 'Pa')
res_u = SaveResults(grid, "u.txt", destinationFolder, 'Displacement', 'm')
res_u.copySettings(sourceFolder, destinationFolder)
# -----------------------------------------------------

# ---------------- INITIAL FIELDS ---------------------
p_old = ScalarField(grid.getNumberOfVertices(), 0.0)
u_old = ScalarField(grid.getNumberOfVertices(), 0.0)
# -----------------------------------------------------

# ------------- CREATE LINEAR SYSTEM ------------------
ls = LinearSystem(2*grid.getNumberOfVertices())
pShift = 1
uShift = (1-pShift)
n = grid.getNumberOfVertices()
# -----------------------------------------------------

# --------------- FLUID FLOW MODEL --------------------
AssemblyDarcyVelocitiesToMatrix(ls, grid, mu, k, pShift)
AssemblyDarcyVelocitiesToVector(ls, grid, mu, k, rho_f, g, pShift)
AssemblyBiotAccumulationToMatrix(ls, grid, timeStep, biot, phi, cs, cf, pShift)
AssemblyBiotAccumulationToVector(ls, grid, timeStep, biot, phi, cs, cf, p_old, pShift)
AssemblyVolumetricStrainToMatrix(ls, grid, timeStep, biot, pShift)
AssemblyVolumetricStrainToVector(ls, grid, timeStep, biot, u_old, pShift)
# -----------------------------------------------------

# -------------- GEOMECHANICAL MODEL ------------------
AssemblyStiffnessMatrix(ls, grid, M, uShift)
AssemblyGravityToVector(ls, grid, rho_s, g, uShift)
AssemblyPorePressureToMatrix(ls, grid, biot, uShift)
u_bottom = 0.0
top_stress = 1.0e3
ls.applyNeumann(n-1 + uShift*n, top_stress)
ls.applyDirichlet(0 + uShift*n, u_bottom)
# -----------------------------------------------------


# --------------------- SOLVE -------------------------
ls.solve()
u_temp, p_temp = ls.splitSolution(n)
print p_temp
p_old.setField(p_temp)
u_old.setField(u_temp)
res_p.saveField(timeHandler.getCurrentTime(), p_old.getField())
res_u.saveField(timeHandler.getCurrentTime(), u_old.getField())
res_p.close()
res_u.close()
# -----------------------------------------------------

print ls.getMatrix()
print ls.getVector()
plotMatrix(ls.getMatrix())




import matplotlib.pyplot as plt
folderName = "results\\"
times = [0]

plt.figure(figsize=(15,6))
plt.subplots_adjust(left=0.05, right=0.95, bottom=0.1, top=0.95)

res_p = ReadResults(folderName + "p.txt")
plt.subplot(1,2,1)
for t in times:
	plt.plot(res_p.getSolutionAtTime(res_p.times[t]), res_p.coord, 'o-')
plt.xlabel("Pressure")
plt.ylabel("Height")
plt.grid(True)

res_u = ReadResults(folderName + "u.txt")
plt.subplot(1,2,2)
for t in times:
	plt.plot(res_u.getSolutionAtTime(res_p.times[t]), res_u.coord, 'o-')
plt.xlabel("Displacement")
plt.ylabel("Height")
plt.grid(True)
plt.show()