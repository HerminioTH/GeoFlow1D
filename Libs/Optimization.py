from UtilitiesLib import computeMedia
import numpy as np

class CellOptimization(object):
	def __init__(self):
		self.cells = []

	def addOptimizer(self, optimizer):
		self.cells.append(optimizer)

	def initialize(self):
		self.size = len(self.cells)
		self.deltas = np.ones(self.size)
		for i in range(self.size):
			self.deltas[i] = self.cells[i].delta

	def optimize(self, ratesOnCells):
		for i in range(self.size):
			self.deltas[i] = self.cells[i].computeDelta(ratesOnCells[i])
		return self.deltas




class BackForwardStep(object):
	def __init__(self, delta_initial, step, factor=10., restart=True):
		self.delta = delta_initial
		self.deltas = [self.delta]
		self.rates = []
		self.initialStep = step
		self.step = self.initialStep
		self.factor = factor
		self.restart = restart

	def computeDelta(self, rate):#, media1, media2):
		self.__saveRate(rate)
		if len(self.rates) < 2:
			return self.delta
		else:
			# if abs(self.rates[-1] - self.rates[-2]) < 1e-2:
			# 	self.step = self.initialStep
			if self.rates[-1] > self.rates[-2]:
				self.__increase()
			else:
				self.__reverse()
				self.__increase()
		return self.delta

	def __increase(self):
		self.delta += self.step
		self.__saveDelta()

	def __reverse(self):
		self.step = -self.step
		self.__reduceStep()

	def __reduceStep(self):
		self.step = self.step/self.factor

	def __saveDelta(self):
		self.deltas.append(self.delta)

	def __saveRate(self, rate):
		if len(self.rates) > 2:
			self.rates.pop(0)
		self.rates.append(rate)



class BackTracking(object):
	def __init__(self, delta_initial, step, factor=10., restart=True):
		self.delta = delta_initial
		self.deltas = [self.delta]
		self.rates = []
		self.initialStep = step
		self.step = self.initialStep
		self.factor = factor
		self.restart = restart

	def computeDelta(self, rate):
		self.__saveRate(rate)
		if len(self.rates) < 2:
			return self.delta
		else:
			if abs(self.rates[-1] - self.rates[-2]) < 1e-2:
				self.step = self.initialStep
			if self.rates[-1] > self.rates[-2]:
				self.__increase()
			else:
				self.__reverse()
				self.__increase()
		return self.delta

	def __increase(self):
		self.delta += self.step
		self.__saveDelta()

	def __reverse(self):
		self.step = -self.step
		self.__reduceStep()

	def __reduceStep(self):
		self.step = self.step/self.factor

	def __saveDelta(self):
		self.deltas.append(self.delta)

	def __saveRate(self, rate):
		if len(self.rates) > 2:
			self.rates.pop(0)
		self.rates.append(rate)