import unittest
from GridLib import *


class Test_Vertex(unittest.TestCase):
	def setUp(self):
		self.v1 = Vertex(11, 0.4)
		self.v1.addToVolume(0.23)
		self.v1.addToVolume(0.31)

	def test_Coordinate(self):
		self.assertEqual(self.v1.getCoordinate(), 0.4)

	def test_AddToVolume(self):
		self.assertEqual(self.v1.getVolume(), 0.54)

	def test_Index(self):
		self.assertEqual(self.v1.getIndex(), 11)

class Test_Face(unittest.TestCase):
	def setUp(self):
		v1 = Vertex(11, 0.4)
		v2 = Vertex(12, 0.5)
		self.face = Face([v1, v2], 2.3)

	def test_NeighborVertices(self):
		self.assertEqual(self.face.getBackwardVertex().getIndex(), 11)
		self.assertEqual(self.face.getForwardVertex().getIndex(), 12)

	def test_Area(self):
		self.assertEqual(self.face.getArea(), 2.3)

	def test_Coordinate(self):
		self.assertEqual(self.face.getCoordinate(), 0.45)

class Test_Element(unittest.TestCase):
	def setUp(self):
		v1 = Vertex(11, 0.4)
		v2 = Vertex(12, 0.52)
		self.e = Element([v1, v2], 3, 2.3)

	def test_Index(self):
		self.assertEqual(self.e.getIndex(), 3)

	def test_Area(self):
		self.assertEqual(self.e.getArea(), 2.3)

	def test_Length(self):
		self.assertEqual(self.e.getLength(), 0.12)

	def test_Volume(self):
		self.assertEqual(self.e.getVolume(), 2.3*0.12)

	def test_Vertices(self):
		self.assertEqual(self.e.getVertices()[0].getIndex(), 11)
		self.assertEqual(self.e.getVertices()[1].getIndex(), 12)

	def test_Face(self):
		self.assertEqual(self.e.getFace().getCoordinate(), (0.4 + 0.52)/2)


class Test_Region(unittest.TestCase):
	def setUp(self):
		v1 = Vertex(11, 0.40)
		v2 = Vertex(12, 0.52)
		v3 = Vertex(17, 0.60)
		e1 = Element([v1, v2], 3, 2.3)
		e2 = Element([v2, v3], 5, 2.3)
		self.region = Region("BODY", 7)
		self.region.addElement(e1)
		self.region.addElement(e2)

	def test_Index(self):
		self.assertEqual(self.region.getIndex(), 7)

	def test_Name(self):
		self.assertEqual(self.region.getName(), "BODY")

	def test_Elements(self):
		elements = self.region.getElements()
		self.assertEqual(elements[0].getIndex(), 3)
		self.assertEqual(elements[1].getIndex(), 5)



class Test_GridData(unittest.TestCase):
	def setUp(self):
		nodesCoord, elemConn = createGridData(10, 11)
		self.gd_1R = GridData()
		self.gd_1R.setElementConnectivity(elemConn)
		self.gd_1R.setNodeCoordinates(nodesCoord)
		self.gd_1R.initialize()

		self.gd_2R = GridData()
		self.gd_2R.setElementConnectivity(elemConn)
		self.gd_2R.setNodeCoordinates(nodesCoord)
		self.gd_2R.setElementsToRegion([0, 1, 2, 3], 'lower_layer')
		self.gd_2R.setElementsToRegion([4, 5, 6, 7, 8, 9], 'upper_layer')
		self.gd_2R.initialize()

	def test_RegionNames(self):
		self.assertEqual(self.gd_1R.regionNames[0], 'None')
		self.assertEqual(self.gd_2R.regionNames[0], 'lower_layer')
		self.assertEqual(self.gd_2R.regionNames[1], 'upper_layer')

	def test_ElementConnectivity(self):
		self.assertListEqual(self.gd_1R.elemConnectivity, [[0,1],[1,2],[2,3],[3,4],[4,5],[5,6],[6,7],[7,8],[8,9],[9,10]])
		self.assertListEqual(self.gd_2R.elemConnectivity, [[0,1],[1,2],[2,3],[3,4],[4,5],[5,6],[6,7],[7,8],[8,9],[9,10]])

	def test_RegionElements(self):
		self.assertListEqual(self.gd_1R.regionElements[0], [0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
		self.assertListEqual(self.gd_2R.regionElements[0], [0, 1, 2, 3])
		self.assertListEqual(self.gd_2R.regionElements[1], [4, 5, 6, 7, 8, 9])

	def test_NodeCoordinates(self):
		self.assertListEqual(list(self.gd_1R.nodeCoordinates), [0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10])
		self.assertListEqual(list(self.gd_2R.nodeCoordinates), [0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10])


class Test_Grid(unittest.TestCase):
	def setUp(self):
		nodesCoord, elemConn = createGridData(15, 11)
		gridData = GridData()
		gridData.setElementConnectivity(elemConn)
		gridData.setNodeCoordinates(nodesCoord)
		gridData.setElementsToRegion([0, 1, 2, 3], 'lower_layer')
		gridData.setElementsToRegion([4, 5, 6, 7, 8, 9], 'upper_layer')
		self.grid = Grid_1D(gridData)

	def test_Numbers(self):
		self.assertEqual(self.grid.getNumberOfVertices(), 11)
		self.assertEqual(self.grid.getNumberOfElements(), 10)
		self.assertEqual(self.grid.getNumberOfRegions(), 2)

	def test_Volumes(self):
		vertices = self.grid.getVertices()
		self.assertEqual(vertices[0].getVolume(), 0.75)
		self.assertEqual(vertices[-1].getVolume(), 0.75)
		for v in vertices[1:-1]:
			self.assertEqual(v.getVolume(), 1.5)

	def test_Regions(self):
		regions = self.grid.getRegions()
		elems_r1 = regions[0].getElements()
		elem_r1_indices = [0, 1, 2, 3]
		[self.assertEqual(elems_r1[i].getIndex(), elem_r1_indices[i]) for i in range(4)]
		elems_r2 = regions[1].getElements()
		elem_r2_indices = [4, 5, 6, 7, 8, 9]
		[self.assertEqual(elems_r2[i].getIndex(), elem_r2_indices[i]) for i in range(4)]

if __name__ == '__main__':
	unittest.main()