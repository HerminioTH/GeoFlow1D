from ResultsHandlerLib import *
import matplotlib.pyplot as plt
# from TerzaghiWithGravity import Solution
from AnalyticalSolution import Solution
from UtilitiesLib import getJsonData

folderName = "results\\"
times = [0, 2, 5, 15, 40, -1]

res_p = ReadResults(folderName + "p.txt")
res_u = ReadResults(folderName + "u.txt")



fluid = getJsonData(folderName + "\\fluid.json")
solid = getJsonData(folderName + "\\solid.json")
H = res_p.coord[-1]
tao = 1e3
g = 0.0
terza = Solution(H, tao, solid, fluid)
np = 50
z_terza = terza.getPositionValues(np)

plt.subplot(1,2,1)
for t in times:
	plt.plot(res_p.getSolutionAtTime(res_p.times[t]), res_p.coord, 'bo')
	plt.plot(terza.getPressureValuesConstTime(res_p.times[t], ny=np), z_terza, 'k-')
plt.grid(True)

plt.subplot(1,2,2)
for t in times:
	plt.plot(res_u.getSolutionAtTime(res_p.times[t]), res_u.coord, 'bo')
	plt.plot(terza.getDisplacementValuesConstTime(res_p.times[t], ny=np), z_terza, 'k-')


plt.grid(True)
plt.show()