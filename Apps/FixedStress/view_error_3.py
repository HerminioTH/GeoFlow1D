
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from UtilitiesLib import *
from ResultsHandlerLib import *
import numpy as np

def plotError(res, ax, times):
	for t in times:
		error = res.getSolutionAtTime(res.times[t])
		ite = range(len(error))
		ax.semilogy(ite, error, 'o-')
	ax.set_xlabel("Iteration number")
	ax.set_ylabel("Error")
	ax.grid(True)

def plotRate(res, ax, name):
	rate = []
	for t in res.times:
		error = res.getSolutionAtTime(t)
		rate.append(computeRate(error))
	# ax.semilogy(res.times, rate, 'o-')
	ax.plot(res.times, rate, '.-', label=name)
	ax.set_xlabel("Time")
	ax.set_ylabel("Rate")
	ax.grid(True)
	ax.legend(loc=0)

def plotIterations(res, ax, name):
	ite = []
	for t in res.times:
		error = res.getSolutionAtTime(t)
		ite.append(len(error))
	# ax.semilogy(res.times, rate, 'o-')
	ax.plot(res.times, ite, '.-', label=name)
	ax.set_xlabel("Time")
	ax.set_ylabel("Number of Iterations")
	ax.grid(True)
	ax.legend(loc=0)

def plotMediaMovel(res, ax, period):
	rate = []
	for t in res.times:
		error = res.getSolutionAtTime(t)
		rate.append(computeRate(error))
	line = []
	time = []
	for i in range(2, len(rate)):
		time.append((res.times[i] + res.times[i-1])/2.)
		media = computeMedia(rate[:i], period)
		line.append(media)
	ax.plot(time, line, '-')#, label="MM %i"%period)







def main():
	# names = ["Case_1\\FIXED_STRAIN\\",
	# 		"Case_3\\FIXED_STRESS_False\\",
	# 		"Case_3\\FIXED_STRESS_True\\"]
	names = ["Case_3\\FIXED_STRESS_False_K\\",
			"Case_3\\FIXED_STRESS_False_M\\",
			"Case_3\\FIXED_STRESS_True_v1\\",
			"Case_3\\FIXED_STRESS_True_v2\\",
			"Case_3\\FIXED_STRESS_True\\"]
	nNames = len(names)

	folderName = "results\\"
	times = [1, 5, 10, 20, -1]
	period = 50

	fig = plt.figure(figsize=(18,10))
	fig.subplots_adjust(left=0.060, right=0.975, top=0.965, bottom=0.08, )
	gs = gridspec.GridSpec(nNames, 3)
	axIte = fig.add_subplot(gs[:,-1])

	for i,name in enumerate(names):
		print i, name
		ax1 = fig.add_subplot(gs[i,0])
		ax2 = fig.add_subplot(gs[i,1])
		res_e = ReadResults(folderName + name + "error.txt")
		plotError(res_e, ax1, times)
		plotRate(res_e, ax2, name)
		plotMediaMovel(res_e, ax2, 3)
		plotMediaMovel(res_e, ax2, 5)
		plotMediaMovel(res_e, ax2, 10)
		plotIterations(res_e, axIte, name)


	# for i,name in enumerate(names):
	# 	print name
	# 	res_e = ReadResults(folderName + name + "error.txt")
	# 	plotError(res_e, axes[i][0], times)
	# 	plotRate(res_e, axes[i][1])
	# 	plotIterations(res_e, axes[0][2])
	# 	plotIterations(res_e, axes[i][2])
	# 	plotMediaMovel(res_e, axes[i][1], 5)
	# 	plotMediaMovel(res_e, axes[i][1], 10)
	# 	plotMediaMovel(res_e, axes[i][1], 20)

	plt.legend(loc=0)
	plt.show()

if __name__ == '__main__':
	main()



