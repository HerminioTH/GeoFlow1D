import sys
import unittest

class IncludeLibsToPath():
	def __init__(self):
		self.directories = []

	def addDirectory(self, dir):
		if type(dir) != str:
			raise "Variable dir has to be a string."
		else:
			self.directories.append(dir)
			if sys.path.count(dir) == 0:
				sys.path.append(dir)

	def __del__(self):
		nd = len(self.directories)
		for i in range(nd):
			sys.path.pop(-1)

class TestIncludeLibsToPath(unittest.TestCase):
	def setUp(self):
		self.pathObj = IncludeLibsToPath()
		self.dirNameA = "TestDirA"
		self.dirNameB = "TestDirB"

	def tearDown(self):
		del self.pathObj
		self.assertNotEqual(sys.path[-2], self.dirNameA)
		self.assertNotEqual(sys.path[-1], self.dirNameB)

	def test_addDirectory(self):
		self.pathObj.addDirectory(self.dirNameA)
		self.pathObj.addDirectory(self.dirNameB)
		self.assertEqual(sys.path[-2], self.dirNameA)
		self.assertEqual(sys.path[-1], self.dirNameB)




if __name__ == "__main__":
	unittest.main()
