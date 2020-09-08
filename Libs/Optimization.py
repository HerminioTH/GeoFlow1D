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
			# if abs(self.rates[-1] - self.rates[-2]) < 1e-3:
			# 	print "------------ RESTART --------------"
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



class ClimbingGoatAlgorithm2(object):
	def __init__(self, delta_initial, step, factor=0.5):
		self.delta = delta_initial
		self.deltas = [self.delta]
		self.rates = []
		self.initialStep = step
		self.step = self.initialStep
		self.factor = factor
		self.direction = 1.
		self.n_reverses = 0

	def computeDelta(self, rate):
		self.anxious_goat_point = False
		self.saveFunctionValue(rate)
		if len(self.rates) < 2:
			return self.delta
			# return self.delta, self.anxious_goat_point
		else:
			self.__peacefulGoat()
		return self.delta
		# return self.delta, self.anxious_goat_point

	def __peacefulGoat(self):
		if self.rates[-1] > self.rates[-2]:
			self.__walk()
			self.n_reverses = 0
			self.__increaseRushedStep()
		else:
			self.__reverse()
			self.__angryGoat()
			self.__walk()

	def __angryGoat(self):
		if self.n_reverses >= 2:
			self.__increaseStep()
			self.anxious_goat_point = True

	def isAngryGoatActive(self):
		return self.anxious_goat_point

	def __walk(self):
		self.delta += self.step*self.direction
		self.__saveDelta()

	def __reverse(self):
		self.direction = -self.direction
		self.__decreaseStep()
		self.n_reverses += 1

	def __decreaseStep(self):
		self.step = self.step/self.factor

	def __increaseRushedStep(self):
		self.step = min(self.initialStep, self.step*1.1)

	def __increasePeacefulStep(self):
		self.step = min(self.initialStep, self.step*(self.factor**1.5))

	def __increaseStep(self):
		self.step = min(self.initialStep, self.step*(self.factor**2.5))

	def __saveDelta(self):
		self.deltas.append(self.delta)

	def saveFunctionValue(self, rate):
		if len(self.rates) > 3:
			self.rates.pop(0)
		self.rates.append(rate)

class ClimbingGoatAlgorithm3(object):
	def __init__(self, delta_initial, step, factor=0.5):
		self.delta = delta_initial
		self.deltas = [self.delta]
		self.rates = []
		self.initialStep = step
		self.step = self.initialStep
		self.factor = factor
		self.direction = 1.
		self.n_reverses = 0
		self.n_peaceful_steps = 0

	def computeDelta(self, rate):
		self.anxious_goat_point = False
		self.saveFunctionValue(rate)
		if len(self.rates) < 2:
			self.__walk()
			return self.delta
		else:
			self.__peacefulGoat()
		return self.delta

	def __peacefulGoat(self):
		if self.rates[-2] > self.rates[-3]:
			self.__walk()
			self.n_reverses = 0
			self.n_peaceful_steps += 1
			self.__increaseRushedStep()
		else:
			self.n_peaceful_steps = 0
			self.__reverse()
			self.__angryGoat()
			self.__walk()

	def __angryGoat(self):
		if self.n_reverses >= 2:
			self.__increaseStep()
			self.anxious_goat_point = True

	def isAngryGoatActive(self):
		return self.anxious_goat_point

	def __walk(self):
		self.delta += self.step*self.direction
		self.__saveDelta()

	def __reverse(self):
		self.direction = -self.direction
		self.__decreaseStep()
		self.n_reverses += 1

	def __decreaseStep(self):
		self.step = self.step/self.factor

	def __increaseRushedStep(self):
		self.step = min(self.initialStep, self.step*1.1)

	def __increasePeacefulStep(self):
		self.step = min(self.initialStep, self.step*(self.factor**1.5))

	def __increaseStep(self):
		self.step = min(self.initialStep, self.step*(self.factor**2.5))

	def __saveDelta(self):
		self.deltas.append(self.delta)

	def saveFunctionValue(self, rate):
		# if len(self.rates) > 3:
		# 	self.rates.pop(0)
		# if len(self.rates) > 2:
		# 	self.rates = []
		self.rates.append(rate)


class ClimbingGoatAlgorithm2D(object):
	def __init__(self, x_0, y_0, step, factor=0.5):
		self.ls_x = ClimbingGoatAlgorithm3(x_0, step, factor)
		self.ls_y = ClimbingGoatAlgorithm3(y_0, step, factor)
		self.pivot = 'x'
		self.x_new = x_0
		self.y_new = y_0
		self.firstCall = False

	def compute(self, f):
		if self.pivot == 'x':
			self.pivot = 'y'
			if self.firstCall:
				self.ls_y.saveFunctionValue(f)
			self.firstCall = True
			self.x_new = self.ls_x.computeDelta(f)
			self.y_new = self.ls_y.delta
		else:
			self.pivot = 'x'
			self.ls_x.saveFunctionValue(f)
			self.x_new = self.ls_x.delta
			self.y_new = self.ls_y.computeDelta(f)
		return self.x_new, self.y_new















class ClimbingGoatAlgorithm(object):
	def __init__(self, delta_initial, step, factor=0.5):
		self.delta = delta_initial
		self.deltas = [self.delta]
		self.rates = []
		self.initialStep = step
		self.step = self.initialStep
		self.factor = factor
		self.direction = 1.
		self.n_reverses = 0

	def computeDelta(self, rate):
		self.anxious_goat_point = False
		self.__saveFunctionValue(rate)
		if len(self.rates) < 2:
			return self.delta
			# return self.delta, self.anxious_goat_point
		else:
			self.__peacefulGoat()
		return self.delta
		# return self.delta, self.anxious_goat_point

	def __peacefulGoat(self):
		if self.rates[-1] > self.rates[-2]:
			self.__walk()
			self.n_reverses = 0
		else:
			self.__reverse()
			self.__angryGoat()
			self.__walk()

	def __angryGoat(self):
		if self.n_reverses >= 2:
			self.__increaseStep()
			self.anxious_goat_point = False

	def isAngryGoatActive(self):
		return self.anxious_goat_point

	def __walk(self):
		self.delta += self.step*self.direction
		self.__saveDelta()

	def __reverse(self):
		self.direction = -self.direction
		self.__decreaseStep()
		self.n_reverses += 1

	def __decreaseStep(self):
		self.step = self.step/self.factor

	def __increaseStep(self):
		self.step = min(self.initialStep, self.step*(self.factor**2.5))

	def __saveDelta(self):
		self.deltas.append(self.delta)

	def __saveFunctionValue(self, rate):
		if len(self.rates) > 3:
			self.rates.pop(0)
		self.rates.append(rate)



class ZigZagSearch(object):
	def __init__(self, x_0, y_0, edgeLength, direction_deg, span_deg, decreaseFactor=1.0, increaseFactor=1.0):
		self.x = [x_0]
		self.y = [y_0]
		self.rates = []
		self.decreaseFactor = decreaseFactor
		self.increaseFactor = increaseFactor
		self.edgeLength_0 = edgeLength
		self.edgeLength = edgeLength
		self.direction = np.radians(direction_deg)
		self.span = np.radians(span_deg)
		self.zigzag = -1.
		self.zigzagStepNumber = 1.

		self.localDeriv = np.array([[-1., 1., 0.], [-1., 0., 1.]])
		self.jacobian = np.zeros((2,2))

		self.firstCall = True
		self.stepNumber = 0

	def compute(self, rate):
		if self.firstCall:
			self.addRate(rate)
			self.takeStep()
			self.firstCall = False
		else:
			self.runRoutine(rate)
		return self.x[-1], self.y[-1]

	def runRoutine(self, rate):
		self.stepNumber += 1
		self.addRate(rate)
		if self.stepNumber == 2:
			self.stepNumber = 0
			self.computeDirection()
			# self.zigzagStepNumber *= 1.
		self.takeStep()

	def takeStep(self):
		self.x.append( self.x[-1] + self.edgeLength*np.cos(self.zigzagStepNumber*self.direction + self.zigzag*self.span) )
		self.y.append( self.y[-1] + self.edgeLength*np.sin(self.zigzagStepNumber*self.direction + self.zigzag*self.span) )
		self.zigzag *= -1.

	def computeDirection(self):
		self.jacobian[0][0] = self.localDeriv[0].dot(self.x[-3:])
		self.jacobian[0][1] = self.localDeriv[0].dot(self.y[-3:])
		self.jacobian[1][0] = self.localDeriv[1].dot(self.x[-3:])
		self.jacobian[1][1] = self.localDeriv[1].dot(self.y[-3:])
		self.globalDeriv = np.linalg.inv(self.jacobian).dot(self.localDeriv)
		s = self.globalDeriv.dot(self.rates[-3:])
		s_norm = np.linalg.norm(s)
		s = s/s_norm

		if s[1] < 0:
			direction_aux = 2*(np.pi - np.arccos(s[0])) + np.arccos(s[0])
		else:
			direction_aux = np.arccos(s[0])

		if abs(direction_aux - self.direction) > np.pi/2.:
			self.edgeLength /= self.decreaseFactor
		else:
			self.edgeLength = min(self.edgeLength_0, self.edgeLength*self.increaseFactor)
		self.direction = direction_aux
		# print direction_aux
		# return s


	def addRate(self, rate):
		self.rates.append(rate)


class ZigZagSearchRandom(object):
	def __init__(self, x_0, y_0, edgeLength, direction_deg, span_deg, decreaseFactor=1.0, increaseFactor=1.0):
		self.x = [x_0]
		self.y = [y_0]
		self.rates = []
		self.decreaseFactor = decreaseFactor
		self.increaseFactor = increaseFactor
		self.edgeLength_0 = edgeLength
		self.edgeLength = edgeLength
		self.direction = np.radians(direction_deg)
		self.span = np.radians(span_deg)
		self.zigzag = -1.
		self.zigzagStepNumber = 1.

		self.localDeriv = np.array([[-1., 1., 0.], [-1., 0., 1.]])
		self.jacobian = np.zeros((2,2))

		self.firstCall = True
		self.stepNumber = 0

		self.changePlane = False
		self.planeChangeCounter = 0

	def compute(self, rate):
		if self.firstCall:
			self.addRate(rate)
			self.takeStep()
			self.firstCall = False
		else:
			self.runRoutine(rate)
		return self.x[-1], self.y[-1]

	def runRoutine(self, rate):
		self.stepNumber += 1
		self.addRate(rate)
		if self.stepNumber == 2:
			self.stepNumber = 0
			self.computeDirection()
		self.takeStep()

	def takeStep(self):
		self.x.append( self.x[-1] + self.edgeLength*np.cos(self.zigzagStepNumber*self.direction + self.zigzag*self.span) )
		self.y.append( self.y[-1] + self.edgeLength*np.sin(self.zigzagStepNumber*self.direction + self.zigzag*self.span) )
		self.zigzag *= -1.

	def computeDirection(self):
		self.jacobian[0][0] = self.localDeriv[0].dot(self.x[-3:])
		self.jacobian[0][1] = self.localDeriv[0].dot(self.y[-3:])
		self.jacobian[1][0] = self.localDeriv[1].dot(self.x[-3:])
		self.jacobian[1][1] = self.localDeriv[1].dot(self.y[-3:])
		self.globalDeriv = np.linalg.inv(self.jacobian).dot(self.localDeriv)
		s = self.globalDeriv.dot(self.rates[-3:])
		s_norm = np.linalg.norm(s)
		s = s/s_norm

		if s[1] < 0:
			direction_aux = 2*(np.pi - np.arccos(s[0])) + np.arccos(s[0])
		else:
			direction_aux = np.arccos(s[0])

		if abs(direction_aux - self.direction) > np.pi/2.:
			self.edgeLength /= self.decreaseFactor
			self.planeChangeCounter += 1
			if self.planeChangeCounter > 2:
				self.changePlane = True
		else:
			self.edgeLength = min(self.edgeLength_0, self.edgeLength*self.increaseFactor)
		self.direction = direction_aux


	def addRate(self, rate):
		self.rates.append(rate)

	def getState(self):
		state = []
		# self.stepNumber = 0
		# self.zigzag = -1.
		state.append(self.edgeLength)
		state.append(self.direction)
		state.append(self.stepNumber)
		state.append(self.zigzag)
		return state

	def loadState(self, state):
		self.edgeLength = state[0]
		self.direction = state[1]
		self.stepNumber = state[2]
		self.zigzag = state[3]



class MultiDirectionalSearch(object):
	def __init__(self, nIntervals, searchAlgorithm, maxNumberOfStepsInPlane=5):
		self.nIntervals = nIntervals
		self.searchAlgorithm = searchAlgorithm
		self.maxNumberOfStepsInPlane = maxNumberOfStepsInPlane

		self.initializeDeltas()
		self.buildPlanes()
		self.initialPlaneConditions()
		self.initializeStates()

	def compute(self, rate):
		if self.shoud_I_change_plane():
			self.saveCurrentPlaneState()
			self.chooseAnotherPlane()
			# self.loadPreviousPlaneState()
			self.adjustDOFsAccordingToNewPlane()
		d0, d1 = self.searchAlgorithm.compute(rate)
		self.appendToDeltas(d0, d1)
		return self.deltas[-1]

	# def shoud_I_change_plane(self):
	# 	if self.searchAlgorithm.changePlane:
	# 		self.searchAlgorithm.changePlane = False
	# 		return True
	# 	else:
	# 		return False

	def shoud_I_change_plane(self):
		self.numberOfStepsInPlane += 1
		if self.numberOfStepsInPlane > self.maxNumberOfStepsInPlane:
			self.numberOfStepsInPlane = 1
			return True
		else:
			return False

	def chooseAnotherPlane(self):
		if self.planeNumber == self.nIntervals-1:
			self.planeNumber = 0
		else:
			self.planeNumber += 1
		self.plane = self.planes[self.planeNumber]

	def saveCurrentPlaneState(self):
		self.states[self.planeNumber] = self.searchAlgorithm.getState()

	def loadPreviousPlaneState(self):
		self.searchAlgorithm.loadState(self.states[self.planeNumber])

	def adjustDOFsAccordingToNewPlane(self):
		self.searchAlgorithm.x[-1] = self.deltas[-1][self.plane[0]]
		self.searchAlgorithm.y[-1] = self.deltas[-1][self.plane[1]]

	def appendToDeltas(self, d0, d1):
		self.aux_deltas = self.deltas[-1].copy()
		self.aux_deltas[self.plane[0]] = d0
		self.aux_deltas[self.plane[1]] = d1
		self.deltas.append(self.aux_deltas)

	def buildPlanes(self):
		self.planes = []
		for i in range(self.nIntervals-1):
			self.planes.append([i, i+1])
		self.planes.append([self.nIntervals-1, 0])

	def initialPlaneConditions(self):
		self.planeNumber = 0
		self.plane = self.planes[0]
		self.numberOfStepsInPlane = 0

	def initializeDeltas(self):
		self.deltas = [[1.0 for i in range(self.nIntervals)]]

	def initializeStates(self):
		self.states = [self.searchAlgorithm.getState() for i in range(self.nIntervals)]
