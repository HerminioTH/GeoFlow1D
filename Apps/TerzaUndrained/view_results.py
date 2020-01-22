from ResultsHandlerLib import *
import matplotlib.pyplot as plt

folderName = "results\\"
times = [0, 2, 5, 15, 40, -1]

res_p = ReadResults(folderName + "p.txt")
plt.subplot(1,2,1)
for t in times:
	plt.plot(res_p.getSolutionAtTime(res_p.times[t]), res_p.coord)
plt.grid(True)

res_u = ReadResults(folderName + "u.txt")
plt.subplot(1,2,2)
for t in times:
	plt.plot(res_u.getSolutionAtTime(res_p.times[t]), res_u.coord)

plt.grid(True)
plt.show()