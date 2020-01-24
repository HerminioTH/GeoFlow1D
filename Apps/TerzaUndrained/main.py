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
nVertices = 7
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
g = -10.0

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

rho = ScalarField(grid.getNumberOfRegions())
for r in grid.getRegions():
	rho.setValue(r, phi.getValue(r)*rho_f + (1 - phi.getValue(r))*rho_s.getValue(r))
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
AssemblyBiotAccumulationToMatrix(ls, grid, timeStep, biot, phi, cs, cf, pShift)
AssemblyVolumetricStrainToMatrix(ls, grid, timeStep, biot, pShift)
# -----------------------------------------------------

# -------------- GEOMECHANICAL MODEL ------------------
AssemblyStiffnessMatrix(ls, grid, M, uShift)
AssemblyPorePressureToMatrix(ls, grid, biot, uShift)
u_bottom = 0.0
top_stress = 1e3
ls.applyDirichletToMatrix(0+uShift*n, u_bottom)
# -----------------------------------------------------



# -------------- TRANSIENT SOLUTION -------------------
timeHandler.advanceTime()
while timeHandler.isFinalTimeReached():
	ls.eraseVector()

	# --------------- FLUID FLOW MODEL --------------------
	AssemblyDarcyVelocitiesToVector(ls, grid, mu, k, rho_f, g, pShift)
	AssemblyBiotAccumulationToVector(ls, grid, timeStep, biot, phi, cs, cf, p_old, pShift)
	AssemblyVolumetricStrainToVector(ls, grid, timeStep, biot, u_old, pShift)
	# -----------------------------------------------------

	# -------------- GEOMECHANICAL MODEL ------------------
	AssemblyGravityToVector(ls, grid, rho, g, uShift)
	ls.applyNeumann(n-1+uShift*n, top_stress)
	ls.applyDirichletToVector(0+uShift*n, u_bottom)
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




import matplotlib.pyplot as plt
from TerzaghiWithGravity import Solution

def getTerza(folderName, H, tao, g):
	fluid = getJsonData(folderName + "\\fluid.json")
	solid = getJsonData(folderName + "\\solid.json")
	return Solution(H, tao, solid, fluid, -g)

folderName = "results\\"
times = [-1]
# times = [0, 10, 50]

res_p = ReadResults(folderName + "p.txt")
res_u = ReadResults(folderName + "u.txt")

terza = getTerza(folderName, L, top_stress, g)
np = 50
z_terza = terza.getPositionValues(np)

plt.figure(figsize=(15,6))
plt.subplots_adjust(left=0.05, right=0.95, bottom=0.1, top=0.95)

plt.subplot(1,2,1)
for t in times:
	plt.plot(res_p.getSolutionAtTime(res_p.times[t]), res_p.coord, 'o')
	plt.plot(terza.getPressureValuesConstTime(0.0, ny=np), z_terza, 'k-')
plt.xlabel("Pressure")
plt.ylabel("Height")
plt.grid(True)

plt.subplot(1,2,2)
for t in times:
	plt.plot(res_u.getSolutionAtTime(res_p.times[t]), res_u.coord, 'o')
	plt.plot(terza.getDisplacementValuesConstTime(0.0, ny=np), z_terza, 'k-')
plt.xlabel("Displacement")
plt.ylabel("Height")
plt.grid(True)
plt.show()


print ls.getVector()

print '\n'
print rho_f*g
print (p_new[-1] - p_new[0])/L
print (terza.getPressureValue(L, 0) - terza.getPressureValue(0, 0))/L

for region in grid.getRegions():
	Q = 1/( cs.getValue(region)*(biot.getValue(region) - phi.getValue(region)) + cf*phi.getValue(region) )
	p_eq = biot.getValue(region)*Q*top_stress/(M.getValue(region) + Q*biot.getValue(region)**2)
print "Equilibrium pressure is: %f Pa"%p_eq