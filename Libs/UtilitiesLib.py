from FieldsLib import *
from GridLib import *
import numpy as np
import matplotlib.pylab as plt
import json

def getJsonData(jsonFile):
    f = open(jsonFile)
    data = json.load(f)
    f.close()
    return data

def computeVolumetricStrain( disp, grid ):
    epsilon = ScalarField( grid.getNumberOfVertices() )
    u = ScalarField( grid.getNumberOfVertices() )
    for i,v in enumerate(grid.getVertices()):
        u.setValue( v, disp[i] )
    for e in grid.getElements():
        f = e.getFace()
        A = f.getArea()
        dx = e.getLength()
        backVertex = f.getBackwardVertex()
        forVertex = f.getForwardVertex()
        operator = [1/(2*dx), 1/(2*dx)]
        operator = [1/2., 1/2.]
        for i,v in enumerate(e.getVertices()):
            value = operator[i]*u.getValue(v)
            epsilon.addValue( backVertex, value/backVertex.getVolume() )
            epsilon.addValue( forVertex, -value/forVertex.getVolume() )
    firstVertex = grid.getVertices()[0]
    lastVertex = grid.getVertices()[-1]
    epsilon.addValue( firstVertex, u.getValue(firstVertex)/firstVertex.getVolume() )
    epsilon.addValue( lastVertex, u.getValue(lastVertex)/lastVertex.getVolume() )
    return epsilon

def updateField( fieldOld, fieldArray, grid ):
    for i,vertex in enumerate(grid.getVertices()):
        fieldOld.setValue(vertex, fieldArray[i])


def plotMatrix( matrix ):
    n = matrix.shape[0]/2
    cmap = plt.cm.rainbow
    fig = plt.figure( 1, (15,10) )
    fig.clf()

    ax = fig.add_subplot( 2, 2, 1 )
    plt.imshow( matrix[0:n,0:n], interpolation='nearest', cmap=cmap )
    plt.colorbar()

    ax = fig.add_subplot( 2, 2, 2 )
    plt.imshow( matrix[0:n,n:], interpolation='nearest', cmap=cmap )
    plt.colorbar()

    ax = fig.add_subplot( 2, 2, 3 )
    plt.imshow( matrix[n:,:n], interpolation='nearest', cmap=cmap )
    plt.colorbar()

    ax = fig.add_subplot( 2, 2, 4 )
    plt.imshow( matrix[n:,n:], interpolation='nearest', cmap=cmap )
    plt.colorbar()
    plt.show()
