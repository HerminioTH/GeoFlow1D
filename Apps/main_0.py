# -*- coding: cp1252 -*-
from FieldsLib import *
from GridLib import *
from LinearSystem import *
# from PropertyParser import Properties
import GeoLib
import numpy as np
import matplotlib.pylab as plt
from PropertyParser import Properties

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
nv = 4
L = 15.
nodesCoord, elemConn = createGridData( L, nv )
gridData = GridData()
gridData.setElementConnectivity( elemConn )
gridData.setNodeCoordinates( nodesCoord )
g = Grid_1D( gridData )
print nodesCoord
## -------------------------------------------------------------------

## -------------------------- SCALAR FIELDS --------------------------
case = "Case_1\\"
rock = Properties( case + "solid.json" )
fluid = Properties( case + "fluid.json" )
time = Properties( case + "NumericalSettings.json" )



# dt = time.getMaterialProperty("Time_Step")

# viscosity = fluid.getMaterialProperty("Viscosity")
# c_f = fluid.getMaterialProperty("Compressibility")

# biot = rock.getMaterialProperty("Biot")
# phi = rock.getMaterialProperty("Porosity")
# c_s = rock.getMaterialProperty("Compressibility")
# poisson = rock.getMaterialProperty("Drained_Poisson")
# shear = rock.getMaterialProperty("Shear_Modulus")

# volumes = ScalarField( g.getNumberOfVertices() )
# g.computeVolumes( volumes )
# permeability = ScalarField( g.getNumberOfElements() )
# lame = ScalarField( g.getNumberOfElements() )

# for e in g.getElements():
#     permeability.setValue( e, rock.getMaterialProperty("Permeability") )
#     lame.setValue( e, 2*shear*poisson/(1-2*poisson) )

# pressure = ScalarField( g.getNumberOfVertices() )
# displacement = ScalarField( g.getNumberOfVertices() )

# equilibriumPressure = pEqui( stress, L, c_f, c_s, phi, poisson, shear, biot )
# print equilibriumPressure,'\n'

# equilibriumDisplacement = uEqui( stress, L, c_f, c_s, phi, poisson, shear, biot )
# print equilibriumDisplacement,'\n'
## -------------------------------------------------------------------