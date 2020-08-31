import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import LightSource
from Optimization import *

def plotShade(Xs, Ts, Fs, ax):
	ls = LightSource(azdeg=15, altdeg=145)
	rgb = ls.shade(Fs, plt.cm.rainbow)
	# cs = ax.contourf(Xs, Ts, Fs, rstride=1, cstride=1, facecolors=rgb, linewidth=0, antialiased=True, shade=True)
	cs = ax.imshow(rgb)

def plotSurface(Xs, Ts, Fs, ax):
	ls = LightSource(azdeg=315, altdeg=45)
	cs = ax.contourf(Xs, Ts, Fs, 50, cmap=cm.rainbow)
	# surf = ax.contour(Xs, Ts, Fs, cmap=cm.rainbow, linewidth=0.1, rstride=1, cstride=1, alpha=1.0)
	# surf = ax.contour(Xs, Ts, Fs, cmap=cm.coolwarm)
	cbar = fig.colorbar(cs)
	cbar.set_label("Rate($\delta$,t)", fontsize=14, **fonts)

def plotLocalMaximum(Xs, Ts, Fs, ax):
	x_max = []
	t_max = []
	# f_max = []
	for i,f_line in enumerate(Fs):
		j = f_line.argmax()
		x_max.append(Xs[i][j])
		t_max.append(Ts[i][j])
		# f_max.append(f_line[j])
	ax.plot(x_max, t_max, '-', color='blue', linewidth=3.5, label='Exhaustive Search')


def plotLine(xs, ts, ax):
	ax.plot(xs, ts, 'k.-', linewidth=1.5, label='The Climbing Goat')

def plotRestartPoints(x, t, ax):
	ax.plot(x,t, 'r.')

def fun(x,t):
	a = -1.8
	b = 5.5
	c = 50.
	# return -(x + a - 1.0*np.sin(-b + t))**2 + 1*t
	return -(x - 2.5*t/(t+5.))**2 + 5*t/(t+5) + 0.1*t

def fun1(x,t):
	sigma = 2 + 30/(t+1)
	mu = 0.032*t + np.sin(0.08*t)
	return np.exp(-0.5*((x - mu)/sigma)**2) / (sigma*np.sqrt(2*np.pi))

def fun2(x,t):
	sigma = 1.5/(0.01*t+1)
	mu = 2 + np.sin(0.1*t)
	return np.exp(-0.5*((x - mu)/sigma)**2) / (sigma*np.sqrt(2*np.pi))

def fun3(x,t):
	sigma = 1.5/(0.01*t+1)
	mu = 1.0 + 1/(1+0.1*t)
	return np.exp(-0.5*((x - mu)/sigma)**2) / (sigma*np.sqrt(2*np.pi))

def fun4(x,t):
	sigma = 1.5/(0.1*np.sin(0.05*t)+1)
	mu = 2 + np.sin(0.1*t)
	return np.exp(-0.5*((x - mu)/sigma)**2) / (sigma*np.sqrt(2*np.pi))

def fun5(x,t):
	sigma = 2.5/(0.1*t+1)
	mu = 2 + np.sin(0.08*t)
	return np.exp(-0.5*((x - mu)/sigma)**2) / (sigma*np.sqrt(2*np.pi))

def fun6(x,t):
	return np.sin(0.6*x)**10 + np.cos(9 + 0.02*t*x) * np.cos(0.9*x)

def fun7(x,t):
	return -(3 - x)**2 - 1*t

if __name__ == "__main__":
	function = fun2
	final_time = 140

	x_initial = 1.0
	x_num = []
	f_num = []
	step = 0.5
	factor = 1.5
	ls = ClimbingGoatAlgorithm2(x_initial, step, factor)
	x_old = x_initial
	restart_x = []
	restart_t = []
	restart_f = []
	anxious_goat_point = 0
	time_num = np.linspace(0, final_time, 140)
	for t in time_num:
		f = function(x_old, t)
		x_num.append(x_old)
		f_num.append(f)
		x_new, anxious_goat_point = ls.computeDelta(f)
		print(anxious_goat_point)
		if anxious_goat_point == 1:
			restart_x.append(x_old)
			restart_f.append(f)
			restart_t.append(t)
		x_old = x_new

	x_min = min(x_num)
	x_max = max(x_num)

	xs = np.linspace(0.0, x_max, 300)
	ts = np.linspace(0, final_time, 300)
	Xs, Ts = np.meshgrid(xs, ts)
	Fs = function(Xs, Ts)



	fig, ax = plt.subplots()
	fig.subplots_adjust(left=0.1, right=1.0, top=0.92, bottom=0.1)
	fonts = {'fontname': 'serif'}
	fontSize = 14

	# ax = fig.gca(projection='2d')
	plotSurface(Xs, Ts, Fs, ax)
	# plotShade(Xs, Ts, Fs, ax)
	plotLocalMaximum(Xs, Ts, Fs, ax)
	plotLine(x_num, time_num, ax)
	plotRestartPoints(restart_x, restart_t, ax)
	ax.set_xlabel("x", fontsize=fontSize, **fonts)
	ax.set_ylabel("Time", fontsize=fontSize, **fonts)
	# ax.set_zlabel("Function")
	# plt.legend(loc=0, shadow=True)
	plt.legend(loc=(0.06, 1.015), fancybox=True, shadow=True, numpoints=1, ncol=5)
	plt.show()

