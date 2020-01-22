from ResultsHandlerLib import *
import pylab as pl

fileName = "p.txt"
folderName = "results\\"
res = ReadResults(folderName + fileName)

times = [res.times[0], res.times[2], res.times[5], res.times[15], res.times[40], res.times[-1]]
for t in times:
	pl.plot(res.coord, res.getSolutionAtTime(t))

pl.grid(True)
pl.show()