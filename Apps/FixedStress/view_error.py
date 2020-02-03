# import matplotlib.pyplot as plt
# from UtilitiesLib import getJsonData
# from ResultsHandlerLib import *

# def computeRate(error):
# 	rate = []
# 	[rate.append(error[i]-error[i-1]) for i in range(1, len(error))]
# 	return rate


# folderName = "results\\FIXED_STRESS\\"
# # folderName = "results\\FIXED_STRAIN\\"
# times = [0, 1, 10, -1]
# # times = [0, 10, 50]

# res_e = ReadResults(folderName + "error.txt")
# print len(res_e.times)


# plt.figure(figsize=(15,6))
# plt.subplots_adjust(left=0.05, right=0.95, bottom=0.1, top=0.95)

# plt.subplot(1,2,1)
# for t in times:
# 	error = res_e.getSolutionAtTime(res_e.times[t])
# 	ite = range(len(error))
# 	plt.semilogy(ite, error, 'o-')
# plt.xlabel("Iteration number")
# plt.ylabel("Error")
# plt.grid(True)
# plt.show()




import matplotlib.pyplot as plt
from UtilitiesLib import getJsonData
from ResultsHandlerLib import *
import numpy as np

# def computeRate(error):
# 	rate = []
# 	[rate.append(error[i-1]-error[i]) for i in range(1, len(error))]
# 	return rate, np.arange(0.5, len(error)-1, 1.0)

def computeRate(error):
	return (np.log10(error[0]) - np.log10(error[-1]))/len(error)

folderName = "results\\FIXED_STRESS_K\\"
# folderName = "results\\FIXED_STRAIN\\"

times = [1, 10, -1]

res_e = ReadResults(folderName + "error.txt")
print len(res_e.times)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16,6))
fig.subplots_adjust(left=0.060, right=0.975, top=0.965, bottom=0.125)
fonts = {'fontname': 'serif'}
fontSize = 14

for t in times:
	error = res_e.getSolutionAtTime(res_e.times[t])
	ite = range(len(error))
	ax1.semilogy(ite, error, 'o-')

rate = []
for t in res_e.times:
	error = res_e.getSolutionAtTime(t)
	rate.append(computeRate(error))
ax2.plot(res_e.times, rate, 'o-')
# ax2.semilogy(res_e.times, rate, 'o-')


ax1.set_xlabel("Iteration number")
ax1.set_ylabel("Error")
ax1.grid(True)

ax2.set_xlabel("Time")
ax2.set_ylabel("Rate")
ax2.grid(True)

plt.show()


