# -*- coding: cp1252 -*-
from FieldsLib import *
from GridLib import *
from LinearSystem import *
from PropertyParser import Properties
import numpy as np
import matplotlib.pylab as plt

def pEqui( stress, height, c_f, c_s, phi, ni, G, alpha ):
    Lambda = 2*G*ni/(1-2*ni)
    Psi = phi*c_f + ( alpha - phi )*c_s
    Beta = Lambda*(1-ni)/ni + alpha*alpha/Psi
    p = -alpha*stress/( Psi*Beta )
    return p

def uEqui( stress, height, c_f, c_s, phi, ni, G, alpha ):
    Psi = phi*c_f + ( alpha - phi )*c_s
    lame = 2*G*ni/(1-2*ni)
    Beta = lame*(1-ni)/ni + alpha*alpha/Psi
    return stress*height/Beta
    
    

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

## -------------------------- GRID --------------------------
def createGridData( L, numberOfNodes ):
    nodesCoord = np.linspace( 0, L, nv )
    elemConn = np.zeros((numberOfNodes-1,2))
    for i in range( numberOfNodes-1 ):
        elemConn[i][0] = i
        elemConn[i][1] = i+1
    return nodesCoord, elemConn

stress = -1e4
nv = 4
L = 15.
nodesCoord, elemConn = createGridData( L, nv )

gridData = GridData()
gridData.setElementConnectivity( elemConn )
gridData.setNodeCoordinates( nodesCoord )

g = Grid_1D( gridData )
## -------------------------------------------------------------------




## -------------------------- SCALAR FIELDS --------------------------
case = "\\Case_0\\"
rock = Properties( case + "Rock_Properties.txt" )
fluid = Properties( case + "Fluid_Properties.txt" )
time = Properties( case + "Numerical_Properties.txt" )

dt = time.getMaterialProperty("Time_Step")

viscosity = fluid.getMaterialProperty("Viscosity")
c_f = fluid.getMaterialProperty("Compressibility")

biot = rock.getMaterialProperty("Biot")
phi = rock.getMaterialProperty("Porosity")
c_s = rock.getMaterialProperty("Compressibility")
poisson = rock.getMaterialProperty("Drained_Poisson")
shear = rock.getMaterialProperty("Shear_Modulus")

volumes = ScalarField( g.getNumberOfVertices() )
g.computeVolumes( volumes )
permeability = ScalarField( g.getNumberOfElements() )
lame = ScalarField( g.getNumberOfElements() )

for e in g.getElements():
    permeability.setValue( e, rock.getMaterialProperty("Permeability") )
    lame.setValue( e, 2*shear*poisson/(1-2*poisson) )

pressure = ScalarField( g.getNumberOfVertices() )
displacement = ScalarField( g.getNumberOfVertices() )

equilibriumPressure = pEqui( stress, L, c_f, c_s, phi, poisson, shear, biot )
print equilibriumPressure,'\n'

equilibriumDisplacement = uEqui( stress, L, c_f, c_s, phi, poisson, shear, biot )
print equilibriumDisplacement,'\n'
## -------------------------------------------------------------------




## -------------------------- LINEAR SYSTEM --------------------------
ls = LinearSystem( 2*g.getNumberOfVertices() )
## -------------------------------------------------------------------



## ------------------------ PRESSURE EQUATION ------------------------
def AssemblyPressureEquation():
    '''Diffusive fluxes'''
    for e in g.getElements():
        dx = e.getLength()
        f = e.getFace()
        A = f.getArea()
        backVertex = f.getBackwardVertex()
        forVertex = f.getForwardVertex()
        bIndex = backVertex.getIndex()
        fIndex = forVertex.getIndex()
        k = permeability.getValue( e )
        diffusiveOperator = [k*A/(viscosity*dx), -k*A/(viscosity*dx)]
    ##    diffusiveOperator = [1, -1]
        localIndex = 0
        for v in e.getVertices():
            flux = diffusiveOperator[localIndex]
            vIndex = v.getIndex()
            ls.addValueToMatrix( bIndex, vIndex, flux )
            ls.addValueToMatrix( fIndex, vIndex, -flux )
            localIndex += 1

    '''Transient terms'''
    for v in g.getVertices():
        psi = phi*c_f + ( biot - phi )*c_s
        V = volumes.getValue( v )
        vIndex = v.getIndex()
        ls.addValueToMatrix( vIndex, vIndex, psi*V/dt )


    '''Volumetric strain terms'''
    for e in g.getElements():
        f = e.getFace()
        A = f.getArea()
        dx = e.getLength()
        backVertex = f.getBackwardVertex()
        forVertex = f.getForwardVertex()
        bIndex = backVertex.getIndex()
        fIndex = forVertex.getIndex()
        operator = [biot*A/(2*dt), biot*A/(2*dt)]
        localIndex = 0
        for v in e.getVertices():
            Vol = volumes.getValue(v)
            flux = operator[localIndex]
            vIndex = v.getIndex() + g.getNumberOfVertices()
            ls.addValueToMatrix( bIndex, vIndex, flux )
            ls.addValueToMatrix( fIndex, vIndex, -flux )
            localIndex += 1
    ls.addValueToMatrix( 0,    nv,     -1*biot*A/(1*dt) )
    ls.addValueToMatrix( nv-1, 2*nv-1, 1*biot*A/(1*dt) )
## -------------------------------------------------------------------




## ----------------------- ELASTICITY EQUATION -----------------------
def AssemblyElasticityEquation():    
    '''Forces'''
    for e in g.getElements():
        dx = e.getLength()
        f = e.getFace()
        A = f.getArea()
        backVertex = f.getBackwardVertex()
        forVertex = f.getForwardVertex()
        bIndex = backVertex.getIndex() + g.getNumberOfVertices()
        fIndex = forVertex.getIndex() + g.getNumberOfVertices()
        value = lame.getValue( e )*(1-poisson)/poisson
        forceOperator = [-value/dx, value/dx]
        localIndex = 0
        for v in e.getVertices():
            flux = forceOperator[localIndex]
            vIndex = v.getIndex() + g.getNumberOfVertices()
            ls.addValueToMatrix( bIndex, vIndex, flux )
            ls.addValueToMatrix( fIndex, vIndex, -flux )
            localIndex += 1

    '''Pore pressure'''
    for e in g.getElements():
        f = e.getFace()
        A = f.getArea()
        backVertex = f.getBackwardVertex()
        forVertex = f.getForwardVertex()
        bIndex = backVertex.getIndex() + g.getNumberOfVertices()
        fIndex = forVertex.getIndex() + g.getNumberOfVertices()
        operator = [biot*A/2., biot*A/2.]
        localIndex = 0
        for v in e.getVertices():
            flux = operator[localIndex]
            vIndex = v.getIndex()
            ls.addValueToMatrix( bIndex, vIndex, flux )
            ls.addValueToMatrix( fIndex, vIndex, -flux )
            localIndex += 1
    ls.addValueToMatrix( 2*nv-1, nv-1, 2*biot*A )
## -------------------------------------------------------------------


AssemblyPressureEquation()
AssemblyElasticityEquation()

## ----------------------- BOUNDARY CONDITIONS -----------------------
##ls.applyDirichlet( 0, 2.0 )
##ls.applyDirichlet( nv-1, 10 )
##ls.applyNeumann( nv-1, 0.0 )

ls.applyDirichlet( nv, 0.0 )
##ls.applyNeumann( nv, -stress )
ls.applyNeumann( 2*nv-1, stress )
## -------------------------------------------------------------------




## ----------------------------- SOLVER ------------------------------
sol = ls.solve()
p, u = splitSolution( sol, g.getNumberOfVertices() )
print "Numerical Solution:"
print p[-1],'\n'
print u[-1]
## -------------------------------------------------------------------

M = ls.getMatrix()
b = ls.getVector()


## ------------------------- POST-PROCESSING -------------------------
x = []
for v in g.getVertices():
    x.append( v.getCoordinate() )

##plt.figure(1)
##plt.plot( p, x )
##plt.xlabel("Pressure")
##plt.ylabel("Position")
##plt.grid(True)
##
##plt.figure(2)
##plt.plot( u, x )
##plt.xlabel("Displacement")
##plt.ylabel("Position")
##plt.grid(True)
##
##plt.show()
## -------------------------------------------------------------------


##plotMatrix( ls.getMatrix() )
