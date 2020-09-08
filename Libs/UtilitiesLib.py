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

def saveDataToJson(fileName, data):
    with open(fileName, "w") as jsonFile:
        json.dump(data, jsonFile, indent=3)

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

def computeMedia(vector, period):
    n = len(vector)
    v = vector[max(n-period, 0):-1]
    if len(v) > 0:
        return sum(v)/len(v)
    else:
        return v[-1]

def computeRate(error):
    return (np.log10(error[0]) - np.log10(max(1e-200, error[-1])))/len(error)
    # return - np.log10(error[-1])/len(error)

# def computeRateOnRegions(error, nIte, grid):


def computeNormL2OnRegions(vector, grid):
    L2 = [[] for i in range(grid.getNumberOfRegions())]
    for region in grid.getRegions():
        soma = 0
        for i,vertex in enumerate(grid.getVerticesFromRegion(region)):
            soma += vertex.getVolume()*vector[i]*vector[i]
        L2[region.getIndex()] = soma**0.5
    return np.array(L2)
