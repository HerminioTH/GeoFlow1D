# -*- coding: cp1252 -*-
from FieldsLib import *
from GridLib import *
from LinearSystem import *
from PropertyParser import Properties
import numpy as np
import matplotlib.pylab as plt
from UtilitiesLib import *

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



## -------------------------- GRID --------------------------
stress = -1e4
nv = 10
L = 15.
nodesCoord, elemConn = createGridData( L, nv )

gridData = GridData()
gridData.setElementConnectivity( elemConn )
gridData.setNodeCoordinates( nodesCoord )

g = Grid_1D( gridData )
## -------------------------------------------------------------------




## -------------------------- SCALAR FIELDS --------------------------
rock = Properties( "Rock_Properties.txt" )
fluid = Properties( "Fluid_Properties.txt" )
time = Properties( "Numerical_Properties.txt" )
timeStep = time.getMaterialProperty("Time_Step")
print timeStep
finalTime = time.getMaterialProperty("Final_Time")
initialTime = 0.0

viscosity = fluid.getMaterialProperty("Viscosity")
c_f = fluid.getMaterialProperty("Compressibility")

biot = rock.getMaterialProperty("Biot")
phi = rock.getMaterialProperty("Porosity")
c_s = rock.getMaterialProperty("Compressibility")
poisson = rock.getMaterialProperty("Drained_Poisson")
shear = rock.getMaterialProperty("Shear_Modulus")

psi = phi*c_f + ( biot - phi )*c_s

permeability = ScalarField( g.getNumberOfElements() )
lame = ScalarField( g.getNumberOfElements() )

for e in g.getElements():
    # permeability.setValue( e, rock.getMaterialProperty("Permeability") )
    lame.setValue( e, 2*shear*poisson/(1-2*poisson) )

for e in g.getElements():
    f = e.getFace()
    centroid = f.getCoordinate()
    if centroid > L/2.:
        permeability.setValue( e, rock.getMaterialProperty("Permeability") )
    else:
        # permeability.setValue( e, rock.getMaterialProperty("Permeability")/1000 )
        permeability.setValue( e, 0.0 )

pressure = ScalarField( g.getNumberOfVertices() )
displacement = ScalarField( g.getNumberOfVertices() )

epsilon_old = ScalarField( g.getNumberOfVertices() )
pressure_old = ScalarField( g.getNumberOfVertices() )

equilibriumPressure = pEqui( stress, L, c_f, c_s, phi, poisson, shear, biot )
equilibriumDisplacement = uEqui( stress, L, c_f, c_s, phi, poisson, shear, biot )
print "Analytical solution:"
print equilibriumPressure
print equilibriumDisplacement,'\n'
## -------------------------------------------------------------------




## -------------------------- LINEAR SYSTEM --------------------------
ls = LinearSystem( 2*g.getNumberOfVertices() )
## -------------------------------------------------------------------



## ------------------------ PRESSURE EQUATION ------------------------
def AssemblyInterpolationMatrix():
    '''Interpolation Function'''
    for e in g.getElements():
        f = e.getFace()
        A = f.getArea()
        dx = e.getLength()
        Lx = dx/3.0
##        Lx = 0.0
        const = Lx*Lx*biot/shear
        operator = [biot*A*const/(dx*timeStep), -biot*A*const/(dx*timeStep)]
        backVertex = f.getBackwardVertex()
        forVertex = f.getForwardVertex()
        bIndex = backVertex.getIndex()
        fIndex = forVertex.getIndex()
        localIndex = 0
        for v in e.getVertices():
            Vol = v.getVolume()
            flux = operator[localIndex]
            vIndex = v.getIndex()
            ls.addValueToMatrix( bIndex, vIndex, flux )
            ls.addValueToMatrix( fIndex, vIndex, -flux )
            localIndex += 1

def AssemblyPressureMatrix():
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
        for i,v in enumerate(e.getVertices()):
            flux = diffusiveOperator[i]
            vIndex = v.getIndex()
            ls.addValueToMatrix( bIndex, vIndex, flux )
            ls.addValueToMatrix( fIndex, vIndex, -flux )

    '''Transient terms'''
    for v in g.getVertices():
        vIndex = v.getIndex()
        psi = phi*c_f + ( biot - phi )*c_s
        ls.addValueToMatrix( vIndex, vIndex, psi*v.getVolume()/timeStep )

    '''Volumetric strain terms'''
    for e in g.getElements():
        f = e.getFace()
        A = f.getArea()
        dx = e.getLength()
        backVertex = f.getBackwardVertex()
        forVertex = f.getForwardVertex()
        bIndex = backVertex.getIndex()
        fIndex = forVertex.getIndex()
        operator = [biot*A/(2*timeStep), biot*A/(2*timeStep)]
        localIndex = 0
        for v in e.getVertices():
            Vol = v.getVolume()
            flux = operator[localIndex]
            vIndex = v.getIndex() + g.getNumberOfVertices()
            ls.addValueToMatrix( bIndex, vIndex, flux )
            ls.addValueToMatrix( fIndex, vIndex, -flux )
            localIndex += 1
    ls.addValueToMatrix( 0,    nv,     -1*biot*A/(1*timeStep) )
    ls.addValueToMatrix( nv-1, 2*nv-1, 1*biot*A/(1*timeStep) )

def AssemblyPressureVector( p_old, e_old ):
    '''Transient terms'''
    for v in g.getVertices():
        psi = phi*c_f + ( biot - phi )*c_s
        V = v.getVolume()
        vIndex = v.getIndex()
        ls.addValueToVector( vIndex, p_old.getValue(v)*psi*V/timeStep )
        ls.addValueToVector( vIndex, biot*V*e_old.getValue(v)/timeStep )
## -------------------------------------------------------------------




## ----------------------- ELASTICITY EQUATION -----------------------
def AssemblyElasticityMatrix():    
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
            flux = -operator[localIndex]
            vIndex = v.getIndex()
            ls.addValueToMatrix( bIndex, vIndex, flux )
            ls.addValueToMatrix( fIndex, vIndex, -flux )
            localIndex += 1
## -------------------------------------------------------------------




solution = TransientSolutionHandler( initialTime, timeStep )




## --------------------------- EQUILIBRIUM ---------------------------
# AssemblyInterpolationMatrix()
AssemblyPressureMatrix()
AssemblyPressureVector( pressure_old, epsilon_old )
AssemblyElasticityMatrix()

## BC Elasticity
ls.applyNeumann( 2*nv-1, stress )
ls.applyDirichlet( nv, 0.0 )

##ls.applyDirichlet( nv-1, 0.0 )

ls.solve()
solution.saveSolution( ls.getSolution() )
p, u = splitSolution( ls.getSolution(), g.getNumberOfVertices() )
print "Numerical Solution:"
print p[-1]
print u[-1],'\n'

epsilon_old = computeVolumetricStrain( u, g )
updateField( pressure_old, p, g )
## -------------------------------------------------------------------




## ---------------------------- TRANSIENT ----------------------------
currentTime = initialTime
while currentTime < timeStep:
    ls.eraseVector()
    AssemblyPressureVector( pressure_old, epsilon_old )

    ## BC Elasticity
    ls.applyNeumann( 2*nv-1, stress )
    ls.applyDirichlet( nv, 0.0 )

    ## BC Pressure
    ls.applyDirichlet( nv-1, 0.0 )

    ls.solve()
    solution.saveSolution( ls.getSolution() )
    p, u = splitSolution( ls.getSolution(), g.getNumberOfVertices() )

    currentTime += timeStep
## -------------------------------------------------------------------



print solution.getPressureOnTime(timeStep)[0]


## ------------------------- POST-PROCESSING -------------------------
x = []
for vertex in g.getVertices():
    x.append( vertex.getCoordinate() )

plt.figure(1)
plt.plot( solution.getPressureOnTime(timeStep), x, 'o-' )
plt.xlabel("Pressure")
plt.ylabel("Position")
plt.grid(True)

##plt.figure(2)
##plt.plot( solution.getDisplacementOnTime(timeStep), x, 'o-' )
##plt.xlabel("Displacement")
##plt.ylabel("Position")
##plt.grid(True)

plt.show()
## -------------------------------------------------------------------


##plotMatrix( ls.getMatrix() )
