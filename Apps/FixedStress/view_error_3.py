
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from UtilitiesLib import *
from ResultsHandlerLib import *
import numpy as np

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

def plotError(res, ax, times):
	for t in times:
		error = res.getSolutionAtTime(res.times[t])
		ite = range(len(error))
		ax.semilogy(ite, error, '.-')
	ax.set_xlabel("Iteration number")
	ax.set_ylabel("Error")
	ax.grid(True)

def plotRate(res, ax, fileName):
	rate = []
	for t in res.times:
		error = res.getSolutionAtTime(t)
		rate.append(computeRate(error))
	pos = fileName.find("STR")
	name = fileName[pos+7:-1]
	# ax.semilogy(res.times, rate, 'o-')
	ax.plot(res.times, rate, '.-', label=name)
	ax.set_xlabel("Time")
	ax.set_ylabel("Rate")
	ax.grid(True)
	ax.legend(loc=0)

def plotIterations(res, ax, fileName):
	ite = []
	for t in res.times:
		error = res.getSolutionAtTime(t)
		ite.append(len(error))
	pos = fileName.find("STR")
	name = fileName[pos+7:-1]
	# ax.loglog(res.times, ite, '.-', label=name)
	ax.semilogy(res.times, ite, '.-', label=name)
	# ax.plot(res.times, ite, '.-', label=name)
	ax.set_xlabel("Time")
	ax.set_ylabel("Number of Iterations")
	ax.grid(True)
	ax.legend(loc=0)

def plotMAXIterations(res, ax, fileName):
	ite = [0]
	for t in res.times:
		error = res.getSolutionAtTime(t)
		ite.append(ite[-1] + len(error))
	pos = fileName.find("STR")
	name = fileName[pos+7:-1]
	# ax.loglog(res.times, ite[1:], '.-', label=name)
	# ax.semilogy(res.times, ite[1:], '.-', label=name)
	ax.plot(res.times, ite[1:], '.-', label=name)
	ax.set_xlabel("Time")
	ax.set_ylabel("Number of Iterations")
	ax.grid(True)
	ax.legend(loc=0)

def plotDeltas(res, ax, fileName):
	deltas = []
	# for t in res.times:
	[deltas.append(res.getSolutionAtTime(t)[0]) for t in res.times]
	pos = fileName.find("STR")
	name = fileName[pos+7:-1]
	# ax.loglog(res.times, ite[1:], '.-', label=name)
	# ax.semilogy(res.times, ite[1:], '.-', label=name)
	ax.plot(res.times, deltas, '.-', label=name)
	ax.set_xlabel("Time")
	ax.set_ylabel("Delta")
	ax.grid(True)
	ax.legend(loc=0)








def main():
	names = [
				"Case_10_SC\\FIXED_STRESS_K\\",
				"Case_10_SC\\FIXED_STRESS_M\\",
				"Case_10_SC\\FIXED_STRESS_D_0.2_2.0_10_20\\",
				"Case_10_SC\\FIXED_STRESS_D_0.2_2.0_10_20_noRestart\\",
				"Case_10_SC\\FIXED_STRESS_D_0.2_2.0_10_20_Restart\\"
				# "Case_8_WC\\FIXED_STRESS_D_0.2_2.0_10_20_noRestart\\",
				# "Case_8_WC\\FIXED_STRESS_D_0.2_2.0_10_20_Restart\\",
				# "Case_9_WC\\FIXED_STRESS_D_0.2_2.0_10_20\\"
				# "Case_8\\FIXED_STRESS_D_0.2_2.0_3_5\\",
				# "Case_8\\FIXED_STRESS_D_0.2_2.0_3_5_noRestart\\",
				# "Case_8\\FIXED_STRESS_D_0.2_2.0_3_5_Restart\\"
				# "Case_5\\FIXED_STRESS_D_0.2_5.0_100_200\\",
				# "Case_5\\FIXED_STRESS_D_0.5_5.0_5_10\\",
				# "Case_5\\FIXED_STRESS_D_0.2_2.0_5_10\\"
			]
	nNames = len(names)

	folderName = "results\\"
	times = [100, 200, -1]

	fig = plt.figure(figsize=(21,13))
	fig.subplots_adjust(left=0.060, right=0.975, top=0.965, bottom=0.08, )
	nCols = 3
	nRows = nNames
	gs = gridspec.GridSpec(nRows, nCols)
	stop = nRows/2
	axIte = fig.add_subplot(gs[stop:,-1])
	axMAX = fig.add_subplot(gs[:stop,-1])
	axRATE = fig.add_subplot(gs[stop:,-2])
	axDELTA = fig.add_subplot(gs[:stop,-2])

	for i,name in enumerate(names):
		print i, name
		ax1 = fig.add_subplot(gs[i,0])
		# ax2 = fig.add_subplot(gs[i,1])
		res_e = ReadResults(folderName + name + "error.txt")
		res_d = ReadResults(folderName + name + "delta.txt")
		plotError(res_e, ax1, times)
		# plotRate(res_e, ax2, name)
		plotRate(res_e, axRATE, name)
		plotDeltas(res_d, axDELTA, name)
		# plotMediaMovel(res_e, ax2, 10)
		# plotMediaMovel(res_e, ax2, 20)
		# plotMediaMovel(res_e, ax2, 10)
		plotIterations(res_e, axIte, name)
		plotMAXIterations(res_e, axMAX, name)

	plt.legend(loc=0)
	plt.show()

if __name__ == '__main__':
	main()



