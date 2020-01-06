import numpy as np

class ScalarField( object ):
    def __init__( self, size, name=None, unity=None ):
        self.__name = name
        self.__unity = unity
        self.__listOfValues = np.zeros( size )

    def setValue( self, entity, value ):
        i = entity.getIndex()
        self.__listOfValues[i] = value

    def addValue( self, entity, value ):
        i = entity.getIndex()
        self.__listOfValues[i] += value

    def getValue( self, entity ):
        i = entity.getIndex()
        return self.__listOfValues[i]

    def getName( self ):
        return self.__name

    def getUnity( self ):
        return self.__unity

    def getValues( self ):
        return self.__listOfValues

    

if __name__ == '__main__':
    from GridLib import *

    elemConn = np.array([[0,1],[1,2],[2,3],[3,4]])
    nodesCoord = np.array([0., 2., 4., 6., 8.])

    gridData = GridData()
    gridData.setElementConnectivity( elemConn )
    gridData.setNodeCoordinates( nodesCoord )

    g = Grid_1D( gridData )

    volumes = ScalarField( g.getNumberOfVertices() )
    width = ScalarField( g.getNumberOfElements() )
