
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from UtilitiesLib import *
from ResultsHandlerLib import *
import numpy as np

def main():
	poissons_str = ["0.0", "0.05", "0.1", "0.15", "0.2", "0.25", "0.3", "0.35", "0.4", "0.45", "0.499"]
	folder = "results\\POISSON\\FIXED_STRESS_D_"

	fig, ax1 = plt.subplots(1, 1, figsize=(6,5))
	fig.subplots_adjust(left=0.1, right=0.965, top=0.965, bottom=0.1)
	fonts = {'fontname': 'serif'}
	fontSize = 14

	poissons_float = []
	deltas = []
	for i,poisson in enumerate(poissons_str):
		poissons_float.append(float(poisson))
		res_d = ReadResults(folder + poisson + "\\delta.txt")
		deltas.append(res_d.getSolutionAtTime(res_d.times[-1])[0])
	ax1.plot(poissons_float, deltas, 'o-', color="black", linewidth=1.5,  markersize=6, markerfacecolor="None", markeredgecolor="black")
	ax1.grid(True)
	ax1.set_ylabel('Delta', fontsize=fontSize, **fonts)
	ax1.set_xlabel("Poisson\'s Ratio", fontsize=fontSize, **fonts)
	ax1.locator_params(axis='y', nbins=20)
	plt.show()

if __name__ == '__main__':
	main()
