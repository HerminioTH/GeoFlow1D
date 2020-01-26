import numpy as np
import scipy.sparse.linalg as spla

class LinearSystem(object):
    def __init__(self, size):
        self.__size = size
        self.__matrix = np.zeros((size, size))
        self.__vector = np.zeros(size)

    def addValueToMatrix(self, row, col, value):
        self.__matrix[row][col] += value

    def setValueToMatrix(self, row, col, value):
        self.__matrix[row][col] = value

    def getMatrixValue(self, row, col):
        return self.__matrix[row][col]

    def addValueToVector(self, row, value):
        self.__vector[row] += value

    def setValueToVector(self, row, value):
        self.__vector[row] = value

    def getVectorValue(self, row):
        return self.__vector[row]

    def getMatrix(self):
        return self.__matrix

    def getVector(self):
        return self.__vector

    def applyDirichlet(self, row, value):
        self.applyDirichletToMatrix(row, value)
        self.applyDirichletToVector(row, value)

    def applyDirichletToMatrix(self, row, value):
        for col in range( self.__size ):
            self.__matrix[row][col] = 0.0
        self.__matrix[row][row] = 1.0

    def applyDirichletToVector(self, row, value):
        self.__vector[row] = value

    def applyNeumann(self, row, value):
        self.__vector[row] += value

    def applyBoundaryConditionsToMatrix(self, grid, boundSettings, shift):
        n = grid.getNumberOfVertices()
        for bName in boundSettings.keys():
            bound = grid.getBoundary(bName)
            bType = boundSettings.get(bName).get("Type")
            bValue = boundSettings.get(bName).get("Value")
            if bType == "Dirichlet":
                self.applyDirichletToMatrix(bound.getVertex().getIndex() + shift*n, bValue)

    def applyBoundaryConditionsToVector(self, grid, boundSettings, shift):
        n = grid.getNumberOfVertices()
        for bName in boundSettings.keys():
            bound = grid.getBoundary(bName)
            bType = boundSettings.get(bName).get("Type")
            bValue = boundSettings.get(bName).get("Value")
            if bType == "Dirichlet":
                self.applyDirichletToVector(bound.getVertex().getIndex() + shift*n, bValue)
            elif bType == "Neumann":
                self.applyNeumann(bound.getVertex().getIndex() + shift*n, bValue)
            else:
                raise Exception("Boundary type %s is not supported."%bType)

    def solve(self):
       # self.__solution, a = spla.bicg( self.__matrix, self.__vector, tol=1e-9, maxiter=1000 )
        self.__solution = np.linalg.solve(self.__matrix, self.__vector)

    def getSolution(self):
        return self.__solution

    def eraseVector(self):
        self.__vector = np.zeros(self.__size)

    def eraseMatrix(self):
        self.__matrix = np.zeros( (self.__size, self.__size) )

    def resetLinearSystem(self):
        self.eraseVector()
        self.eraseMatrix()

    def splitSolution(self, n):
        return self.__solution[:n], self.__solution[n:]


class TransientSolutionHandler( object ):
    def __init__(self, initialTime, timeStep):
        self.__initialTime = initialTime
        self.__timeStep = timeStep
        self.__time = initialTime - timeStep
        self.__solution = []
        self.__timeList = []

    def saveSolution(self, solution):
        self.__time += self.__timeStep
        self.__timeList.append(self.__time)
        self.__solution.append(solution)

    def getSolutionAtTime(self, time):
        try:
            i = self.__timeList.index(time)
            return self.__solution[i]
        except:
            print "Time %f do not exist."%time

    def getPressureAtTime(self, time):
        try:
            i = self.__timeList.index(time)
            try:    n = self.__solution[i].size/2
            except: n = len(self.__solution[i])/2
            return self.__solution[i][:n]
        except:
            print "Time %f do not exist."%time

    def getDisplacementAtTime(self, time):
        try:
            i = self.__timeList.index(time)
            try:    n = self.__solution[i].size/2
            except: n = len( self.__solution[i] )/2
            return self.__solution[i][n:]
        except:
            print "Time %f do not exist."%time

    def getSolutionOnIndex(self, index=0):
        return self.__solution[index]

    def getPressureOnIndex(self, index=0):
        try:    n = self.__solution[index].size/2
        except: n = len( self.__solution[index] )/2
        return self.__solution[index][:n]

    def getDisplacementOnIndex(self, index=0):
        try:    n = self.__solution[index].size/2
        except: n = len( self.__solution[index] )/2
        return self.__solution[index][n:]

    def getSolution(self):
        return self.__solution

    def getTimeList(self):
        return self.__timeList


def splitSolution(sol, n):
    return sol[:n],sol[n:]


if __name__ == '__main__':
    initialTime = 5.0
    timeStep = 0.2
    s = TransientSolutionHandler(initialTime, timeStep)
    s.saveSolution([1,2,3,4])
    s.saveSolution([5,6,7,8])
    s.getSolutionOnTime(9)
    print s.getSolutionOnTime(5.2)
