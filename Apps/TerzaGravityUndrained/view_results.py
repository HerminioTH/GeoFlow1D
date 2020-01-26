import matplotlib.pyplot as plt
from UtilitiesLib import getJsonData
from TerzaghiWithGravity import Solution
from ResultsHandlerLib import *

def getTerza(folderName, H, tao, g):
	fluid = getJsonData(folderName + "\\fluid.json")
	solid = getJsonData(folderName + "\\solid.json")
	return Solution(H, tao, solid, fluid, -g)

folderName = "results\\"
times = [0, 1, -1]
# times = [0, 10, 50]

res_p = ReadResults(folderName + "p.txt")
res_u = ReadResults(folderName + "u.txt")
L = res_p.coord[-1]

ic = getJsonData(folderName + "IC.json")
g = ic.get("Gravity")
bound_u = getJsonData(folderName + "BC_u.json")
top_stress = bound_u.get("TOP").get("Value")

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