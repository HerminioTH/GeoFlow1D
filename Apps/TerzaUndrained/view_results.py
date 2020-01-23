from ResultsHandlerLib import *
import matplotlib.pyplot as plt
from TerzaghiWithGravity import Solution
from UtilitiesLib import getJsonData

folderName = "results\\"
times = [0]

res_p = ReadResults(folderName + "p.txt")
H = res_p.coord[-1]
fluid = getJsonData(folderName + "\\fluid.json")
solid = getJsonData(folderName + "\\solid.json")
tao = 1e3
g = 0.0
terza = Solution(H, tao, solid, fluid, g)
z_terza = terza.getPositionValues(50)
p_terza = terza.getPressureValue(z_terza, 0.0)
u_terza = terza.getDisplacementValue(z_terza, 0.0)

print p_terza
print res_p.getSolutionAtTime(res_p.times[0])

plt.subplot(1,2,1)
for t in times:
	plt.plot(res_p.getSolutionAtTime(res_p.times[t]), res_p.coord, 'o')
	plt.plot(p_terza, z_terza, '-')
plt.grid(True)

res_u = ReadResults(folderName + "u.txt")
plt.subplot(1,2,2)
for t in times:
	plt.plot(res_u.getSolutionAtTime(res_p.times[t]), res_u.coord, 'o')
	plt.plot(u_terza, z_terza, '-')


plt.grid(True)
plt.show()