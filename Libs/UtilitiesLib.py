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


def computeNormL2(vector, grid):
    soma = 0
    for i,vertex in enumerate(grid.getVertices()):
        soma += vertex.getVolume()*vector[i]*vector[i]
    return soma**0.5

def computeNormInf(vector):
    return vector.max()