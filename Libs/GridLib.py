import numpy as np

class GridData( object ):
    def __init__( self ):
        self.elemConnectivity = []
        self.nodeCoordinates = []

    def setElementConnectivity( self, elemConnectivity ):
        self.elemConnectivity = elemConnectivity

    def setNodeCoordinates( self, nodeCoordinates ):
        self.nodeCoordinates = nodeCoordinates

class Vertex( object ):
    def __init__( self, index, coord ):
        self.__globalIndex = index
        self.__x = coord
        self.__volume = 0.0

    def getIndex( self ):
        return self.__globalIndex
    
    def getCoordinate( self ):
        return self.__x

    def addToVolume( self, value ):
        self.__volume += value

    def getVolume( self ):
        return self.__volume

    

class Face( object ):
    def __init__( self, vertices, height=1.0 ):
        self.__bVertex = vertices[0]
        self.__fVertex = vertices[1]
        self.__faceCoord = ( vertices[0].getCoordinate() + vertices[1].getCoordinate() ) / 2.0
        self.__area = height

    def getBackwardVertex( self ):
        return self.__bVertex

    def getForwardVertex( self ):
        return self.__fVertex

    def getCoordinate( self ):
        return self.__faceCoord

    def getArea( self ):
        return self.__area
    

class Element( object ):
    def __init__( self, vertices, index, height=1.0 ):
        self.__globalIndex = index
        self.__height = height
        self.__vertices = vertices
        self.__buildFace()
        self.__elementLength = self.__vertices[1].getCoordinate() - self.__vertices[0].getCoordinate()

    def __buildFace( self ):
        self.__face = Face( self.__vertices )

    def getIndex( self ):
        return self.__globalIndex

    def getVertices( self ):
        return self.__vertices

    def getFace( self ):
        return self.__face

    def getLength( self ):
        return self.__elementLength

    def getHeight( self ):
        return self.__height


class Grid_1D( object ):
    def __init__( self, gridData ):
        self.__buildVertices( gridData )
        self.__buildElements( gridData )
        self.__computeVolumes()

    def __buildVertices( self, gridData ):
        self.__vertices = []
        iVertices = []        
        for iElem in gridData.elemConnectivity:
            for iVertex in iElem:
                if iVertices.count( iVertex ) == 0:
                    iVertices.append( iVertex )
                    self.__vertices.append( Vertex( iVertex, gridData.nodeCoordinates[iVertex] ) )
        self.__numberOfVertices = len( self.__vertices )
                

    def __buildElements( self, gridData ):
        self.__elements = []
        self.__nElements = 0
        for i, iElem in enumerate(gridData.elemConnectivity):
            self.__elements.append( Element( [ self.__vertices[iElem[0]], self.__vertices[iElem[1]] ], i ) )
            self.__nElements += 1
            

    def getVertices( self ):
        return self.__vertices
            
    def getElements( self ):
        return self.__elements

    def getNumberOfElements( self ):
        return self.__nElements

    def getNumberOfVertices( self ):
        return self.__numberOfVertices

    def __computeVolumes( self ):
        for e in self.__elements:
            subVol = e.getLength()/2.
            for v in e.getVertices():
                v.addToVolume( subVol )


def createGridData( L, numberOfNodes ):
    nodesCoord = np.linspace( 0, L, numberOfNodes )
    elemConn = np.zeros((numberOfNodes-1,2),dtype=int)
    for i in range( numberOfNodes-1 ):
        elemConn[i][0] = i
        elemConn[i][1] = i+1
    return nodesCoord, elemConn


if __name__ == '__main__':
    from FieldsLib import *

    L = 10.
    nVertices = 3
    nodesCoord, elemConn = createGridData( L, nVertices )

    gridData = GridData()
    gridData.setElementConnectivity( elemConn )
    gridData.setNodeCoordinates( nodesCoord )

    v1 = Vertex( 0, 1.5 )
    v2 = Vertex( 1, 1.7 )
    e1 = Element( [v1,v2], 12 )

    g = Grid_1D( gridData )

    
