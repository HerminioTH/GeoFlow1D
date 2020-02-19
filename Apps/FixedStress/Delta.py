from UtilitiesLib import computeMedia

class Delta(object):
	def __init__(self, delta_initial, step, factor=10., restart=True):
		self.delta = delta_initial
		self.deltas = [self.delta]
		self.rates = []
		self.initialStep = step
		self.step = self.initialStep
		self.factor = factor
		self.restart = restart
		self.counter = 0
		self.reduceCounter = 0

	def computeDelta(self, rate, media1, media2):
		self.__saveRate(rate)
		if self.counter < 2:
			self.counter += 1
			return self.delta
		else:
			if self.rates[-1] > self.rates[-2]:
				self.__increase()
			else:
				if len(self.deltas) > 1:
					delta_medio = computeMedia(self.deltas, 10)
				else:
					delta_medio = 2*self.delta
				if self.restart:
					if rate<media1 and rate<media2 and abs(self.delta - delta_medio)<1e-6:
						self.step = self.initialStep
						print "----------------------------------"
						print "----------- RESTART --------------"
						print "-----------" + str(self.step) + "--------------"
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
		self.rates.append(rate)