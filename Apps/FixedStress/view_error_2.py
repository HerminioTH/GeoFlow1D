
import matplotlib.pyplot as plt
from UtilitiesLib import *
from ResultsHandlerLib import *
import numpy as np

# def computeRate(error):
# 	rate = []
# 	[rate.append(error[i-1]-error[i]) for i in range(1, len(error))]
# 	return rate, np.arange(0.5, len(error)-1, 1.0)

def plotError(res, ax, times):
	for t in times:
		error = res.getSolutionAtTime(res.times[t])
		ite = range(len(error))
		ax.semilogy(ite, error, 'o-')
	ax.set_xlabel("Iteration number")
	ax.set_ylabel("Error")
	ax.grid(True)

def plotRate(res, ax):
	rate = []
	for t in res.times:
		error = res.getSolutionAtTime(t)
		rate.append(computeRate(error))
	# ax.semilogy(res.times, rate, 'o-')
	ax.plot(res.times, rate, '.')
	ax.set_xlabel("Time")
	ax.set_ylabel("Rate")
	ax.grid(True)

def plotMediaMovel(res, ax, period):
	rate = []
	for t in res.times:
		error = res.getSolutionAtTime(t)
		rate.append(computeRate(error))
	line = []
	time = []
	print len(rate)
	for i in range(2, len(rate)):
		time.append((res.times[i] + res.times[i-1])/2.)
		media = computeMedia(rate[:i], period)
		line.append(media)
	ax.plot(time, line, '-', label="MM %i"%period)







def main():
	names = ["FIXED_STRAIN\\",
			"FIXED_STRESS_K\\",
			"FIXED_STRESS_M\\"]
	nNames = len(names)

	folderName = "results\\"
	times = [1, 10, 30, -1]
	period = 50

	fig, axes = plt.subplots(nNames, 2, figsize=(20,20))
	fig.subplots_adjust(left=0.060, right=0.975, top=0.965, bottom=0.125)
	fonts = {'fontname': 'serif'}
	fontSize = 14

	for i,name in enumerate(names):
		res_e = ReadResults(folderName + name + "error.txt")
		plotError(res_e, axes[i][0], times)
		plotRate(res_e, axes[i][1])
		plotMediaMovel(res_e, axes[i][1], 5)
		plotMediaMovel(res_e, axes[i][1], 10)
		plotMediaMovel(res_e, axes[i][1], 20)

	plt.legend(loc=0)
	plt.show()

if __name__ == '__main__':
	main()



