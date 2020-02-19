from UtilitiesLib import computeMedia

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