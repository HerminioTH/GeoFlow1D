import unittest
from GridLib import *
from FlowLib import *
from FieldsLib import *
from LinearSystem import *

class Test_Flow(unittest.TestCase):
	def setUp(self):
		L = 6.
		n = 3
		self.dx = L/(n-1.)
		A = 1.0
		nodesCoord, elemConn = createGridData(L, n)
		gridData = GridData()
		gridData.setElementConnectivity(elemConn)
		gridData.setNodeCoordinates(nodesCoord)
		self.grid = Grid_1D(gridData)

		k = 3.0
		self.permeability = ScalarField(self.grid.getNumberOfRegions())
		self.permeability.setValue(self.grid.getRegions()[0], k)
		self.viscosity = 1e-3
		self.density = 1000.
		self.gravity = -10.

		self.D = k*A/self.viscosity
		self.Dx = self.D/self.dx

		self.ls = LinearSystem(self.grid.getNumberOfVertices())

	def test_AssemblyDarcyVelocities(self):
		AssemblyMassDarcyVelocities(self.ls, self.grid, self.viscosity, self.permeability, self.density, self.gravity, pShift=0)
		self.assertEqual(self.ls.getMatrixValue(0,0), self.Dx)
		self.assertEqual(self.ls.getMatrixValue(0,1), -self.Dx)
		self.assertEqual(self.ls.getMatrixValue(1,0), -self.Dx)
		self.assertEqual(self.ls.getMatrixValue(1,1), 2*self.Dx)
		self.assertEqual(self.ls.getMatrixValue(1,2), -self.Dx)
		self.assertEqual(self.ls.getMatrixValue(2,1), -self.Dx)
		self.assertEqual(self.ls.getMatrixValue(2,2), self.Dx)

		self.assertEqual(self.ls.getVectorValue(0), -self.D*self.density*self.gravity)
		self.assertEqual(self.ls.getVectorValue(2), self.D*self.density*self.gravity)







if __name__ == '__main__':
	unittest.main()