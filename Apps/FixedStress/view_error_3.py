
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from UtilitiesLib import *
from ResultsHandlerLib import *
import numpy as np

fonts = {'fontname': 'serif'}
fontSize = 14

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

def plotError(res, ax, times, fileName):
	maxIte = 0
	for t in times:
		error = res.getSolutionAtTime(res.times[t])
		ite = range(1,len(error)+1)
		if ite[-1] > maxIte:
			maxIte = ite[-1]
		ax.semilogy(ite, error, '.-')

	pos = fileName.find("STR")
	name = fileName[pos+7:-1]
	ax.semilogy(ite, error, '.-', label=name)
	ax.set_xlabel("Iteration number", fontsize=fontSize, **fonts)
	ax.set_ylabel("Error", fontsize=fontSize, **fonts)
	ax.locator_params(axis='x', nbins=10)
	ax.set_xlim(1, maxIte+1)
	ax.grid(True)
	ax.legend(loc=0)

def plotRate(res, ax, fileName):
	rate = []
	for t in res.times:
		error = res.getSolutionAtTime(t)
		rate.append(computeRate(error))
	pos = fileName.find("STR")
	name = fileName[pos+7:-1]
	# ax.semilogy(res.times, rate, 'o-')
	ax.plot(res.times, rate, '.-', label=name)
	ax.set_xlabel("Time [s]", fontsize=fontSize, **fonts)
	ax.set_ylabel("Rate", fontsize=fontSize, **fonts)
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
	ax.set_xlabel("Time [s]", fontsize=fontSize, **fonts)
	ax.set_ylabel("Number of Iterations", fontsize=fontSize, **fonts)
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
	ax.set_xlabel("Time [s]", fontsize=fontSize, **fonts)
	ax.set_ylabel("Number of Iterations", fontsize=fontSize, **fonts)
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
	ax.set_xlabel("Time [s]", fontsize=fontSize, **fonts)
	ax.set_ylabel("Delta", fontsize=fontSize, **fonts)
	ax.grid(True)
	ax.legend(loc=0)








def main():
	names = [
				"Case_10_SC\\FIXED_STRESS_K\\",
				"Case_10_SC\\FIXED_STRESS_M\\",
				"Case_10_SC\\FIXED_STRESS_D_noRestart\\",
				"Case_10_SC\\FIXED_STRESS_D_Restart\\"
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
	times = [0, 1, 10, 50, 100, 200, -100, -50, -1]

	fig = plt.figure(figsize=(21,12))
	fig.subplots_adjust(left=0.060, right=0.975, top=0.965, bottom=0.08, hspace=0.3)
	nCols = 3
	nRows = nNames
	gs = gridspec.GridSpec(nRows, nCols)
	stop = nRows/2
	axIte = fig.add_subplot(gs[stop:,-1])
	axMAX = fig.add_subplot(gs[:stop,-1])
	axRATE = fig.add_subplot(gs[stop:,-2])
	axDELTA = fig.add_subplot(gs[:stop,-2])

	for i,name in enumerate(names):
		ax1 = fig.add_subplot(gs[i,0])
		# ax2 = fig.add_subplot(gs[i,1])
		res_e = ReadResults(folderName + name + "error.txt")
		res_d = ReadResults(folderName + name + "delta.txt")

		print i, name, len(res_e.times)
		plotError(res_e, ax1, times, name)
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



