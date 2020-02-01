"""
Created on Thu Nov 7 14:08:02 2019

@author: Carlos
         Adapted from Herminio's code.
"""
import numpy as np
import pylab as pl
import math

# CLASS DEFINITION ==============================================================================

class Solution( object ):
    def __init__( self, height, tao_0, rock, fluid, gravity ):
        self.height = height;
        self.tao_0 = tao_0;
        self.g = gravity;

        fluidName = fluid.keys()[0]
        self.mi = fluid.get(fluidName).get("Viscosity").get("value")
        self.c_f = fluid.get(fluidName).get("Compressibility").get("value")
        self.rho_f = fluid.get(fluidName).get("Density").get("value")

        solidName = rock.keys()[0]
        self.c_s = rock.get(solidName).get("Compressibility").get("value")
        self.permeability = rock.get(solidName).get("Permeability").get("value")
        self.phi = rock.get(solidName).get("Porosity").get("value")
        self.ni = rock.get(solidName).get("PoissonsRatio").get("value")
        self.G = rock.get(solidName).get("ShearModulus").get("value")
        self.rho_s = rock.get(solidName).get("Density").get("value")

        self.rho = self.rho_f * self.phi + ( 1 - self.phi) * self.rho_s
        self.lame1st = 2 * self.G * self.ni / ( 1 - 2 * self.ni )
        self.M = 2 * self.G + self.lame1st

        # self.mi = fluid.fromMaterialGetProperty( 0, 'Viscosity' );
        # self.c_f = fluid.fromMaterialGetProperty( 0, 'Compressibility' );
        # self.rho_f = fluid.fromMaterialGetProperty( 0, 'Density' );

        # self.c_s = rock.fromMaterialGetProperty( 0, 'Compressibility' );
        # self.permeability = rock.fromMaterialGetProperty( 0, 'Permeability' );
        # self.phi = rock.fromMaterialGetProperty( 0, 'Porosity' );
        # self.ni = rock.fromMaterialGetProperty( 0, 'PoissonsRatio' );
        # # self.ni_u = rock.fromMaterialGetProperty( 0, 'UndrainedPoissonsRatio' )
        # self.G = rock.fromMaterialGetProperty( 0, 'ShearModulus' );
        # self.rho_s = rock.fromMaterialGetProperty( 0, 'Density' );
        # self.rho = self.rho_f * self.phi + ( 1 - self.phi) * self.rho_s;
        # self.lame1st = 2 * self.G * self.ni / ( 1 - 2 * self.ni );
        # self.M = 2 * self.G + self.lame1st;

        self._calculate_K_s();
        self._calculate_K_f();
        self._calculate_K_phi();
        self._calculate_K();

        self._calculate_alpha();
        self._calculate_Q();

        self._calculate_H();
        self._calculate_R();
        self._calculate_K_p();
        self._calculate_B();
        self._calculate_K_u();
        self._calculate_ni_u();
        self._calculate_K_ni_u();
        self._calculate_gama();
        self._calculate_c();

        self._calculate_K_v_u();
        self._calculate_c_m();

    # Internal functions --------------------------------------------------------------------

    def _calculate_K_s( self ):
        if self.c_s == 0.0:
            self.K_s = 1.0e+100
        else:
            self.K_s = 1.0 / self.c_s;

    def _calculate_K_f( self ):
        if self.c_f == 0.0:
            self.K_f = 1.0e+100
        else:
            self.K_f = 1.0 / self.c_f;

    def _calculate_K_phi( self ):
        self.K_phi = self.K_s

    def _calculate_K( self ):
        self.K = 2*self.G*( 1 + self.ni ) / ( 3 - 6*self.ni )

    def _calculate_alpha(self):
        self.alpha = 1 - self.c_s*self.K

    def _calculate_H( self ):
        self.H = 1.0 / ( ( 1.0 / self.K ) - ( 1.0 / self.K_s ) );

    def _calculate_R( self ):
        self.R = 1.0 / ( ( 1.0 / self.H ) + self.phi * ( ( 1.0 / self.K_f ) - ( 1.0 / self.K_phi ) ) );

    def _calculate_K_p( self ):
        self.K_p = self.phi * self.K / self.alpha;

    def _calculate_B( self ):
        self.B = self.R / self.H;

    def _calculate_K_u( self ):
        if self.alpha * self.B != 1.0:
            self.K_u = self.K / ( 1.0 - self.alpha * self.B );
        else:
            self.K_u = 1.0e+100

    def _calculate_ni_u( self ):
        self.ni_u = ( ( 3.0 * self.ni + self.alpha * self.B * ( 1.0 - 2.0 * self.ni ) ) /
                          ( 3.0 - self.alpha * self.B * ( 1.0 - 2.0 * self.ni ) ) );

    def _calculate_K_ni_u( self ):
        self.K_ni_u = ( 3.0 * self.K_u * ( 1.0 - self.ni_u ) ) / ( 1.0 + self.ni_u );

    def _calculate_gama( self ):
        self.gama = ( self.B * ( 1.0 + self.ni_u ) ) / ( 3.0 * ( 1.0 - self.ni_u ) );

    def _calculate_c( self ):
        self.c = ( ( 2.0 * self.permeability * self.G * ( 1.0 - self.ni ) * ( self.ni_u - self.ni ) ) /
                    ( self.mi * ( self.alpha ** 2.0 ) * ( 1.0 - self.ni_u ) * ( ( 1.0 - 2.0 * self.ni ) ** 2.0 ) ) );

    def _calculate_K_v_u( self ):
        self.K_v_u = 3*self.K_u*( 1 - self.ni_u ) / ( 1 + self.ni_u )

    def _calculate_c_m( self ):
        self.c_m = ( self.alpha * ( 1.0 + self.ni ) ) / ( 3.0 * self.K * ( 1.0 - self.ni ) );

    def _calculate_ni( self ):
        self.ni = ( 3.0 * self.K - 2.0 * self.G ) / ( 2.0 * ( 3.0 * self.K + self.G ) );

    def _calculate_alpha( self ):
        self.alpha = 1.0 - ( self.K / self.K_s );

    def _calculate_Q( self ):
        self.Q = 1 / ( self.phi * self.c_f + ( self.alpha - self.phi ) * self.c_s );

    def _calculate_p_0( self, yPosition ):
        Mu = self.M + self.alpha*self.alpha*self.Q
        p_0 = self.alpha*self.Q*(self.tao_0 + 0.5*self.rho*self.g*self.height) / Mu \
              - self.rho_f*self.g*(yPosition - 0.5*self.height)
        return p_0

    def _calculate_v_0( self, yPosition ):
        Mu = self.M + self.alpha*self.alpha*self.Q
        v_0 = (self.rho - self.alpha*self.rho_f)*self.g*yPosition*(yPosition - self.height) / (2*self.M) \
              - (self.tao_0 + 0.5*self.rho*self.g*self.height)*yPosition / Mu
        return v_0

    # def getPositionValues( self, n = 200, axisName = None ):
    #     return self.getPositionValues(n)

    # def getPositionValues( self, ny = 200 ):
    #     dy = self.height / ( ny - 1.0 );
    #     positionValues = [];
    #     for i in range( 0, ny ):
    #         positionValues.append( i * dy );
    #     return np.array(positionValues)

    def getPositionValues(self, n=200):
        return np.linspace(0, self.height, n)

    def getPressureValue( self, yPosition, time, numberOfSummationTerms = 200 ):
        position = self.height - yPosition;
        if time == 0.0:
            pressureValue = self._calculate_p_0(yPosition);
            return pressureValue
        else:
            summationResult = 0
            for j in range(0, numberOfSummationTerms):
                term_1 = 1.0 / (2.0*j + 1.0)
                term_2 = math.exp( -(time*self.c*(math.pi**2)*(2.0*j + 1.0)**2.0 ) / (4.0*self.height**2.0) )
                term_3 = math.sin( math.pi*position*(2.0*j + 1) / (2.0*self.height) )
                summationResult += term_1*term_2*term_3
            barP0 = self._calculate_p_0(self.height)
            pressureValue = self.rho_f*self.g*position + 4.0*barP0*summationResult / math.pi
            return pressureValue

    # def getPressureValue( self, yPosition, time, numberOfSummationTerms = 200 ):
    #     position = self.height - yPosition;
    #     if time == 0.0:
    #         pressureValue = self._calculate_p_0(yPosition);
    #         return pressureValue
    #     else:
    #         summationResult = 0
    #         for j in range(0, numberOfSummationTerms):
    #             term_1 = 1.0 / (2.0*j + 1.0)
    #             term_2 = math.exp( -(time*self.c*(math.pi**2)*(2.0*j + 1.0)**2.0 ) / (4.0*self.height**2.0) )
    #             term_3 = math.cos( math.pi*position*(2.0*j + 1) / (2.0*self.height) )
    #             summationResult += term_1*term_2*term_3
    #         barP0 = self._calculate_p_0(self.height)
    #         pressureValue = self.rho_f*self.g*position + 4.0*barP0*summationResult / math.pi
    #         return pressureValue

    def getDisplacementValue( self, yPosition, time, numberOfSummationTerms = 200 ):
        position = self.height - yPosition;
        if time == 0.0:
            displacementValue = self._calculate_v_0( yPosition );
            return displacementValue
        else:
            summationResult = 0;
            for j in range( 0, numberOfSummationTerms ):
                term_1 = 1.0 / ( ( 2.0 * j + 1.0 ) ** 2 )
                term_2 = ( math.exp( - ( ( self.c * time * ( math.pi ** 2.0 ) * ( ( 2.0 * j + 1.0 ) ** 2.0 ) ) /
                                         ( 4.0 * ( self.height ** 2.0 ) ) ) ) )
                term_3 = math.cos( ( math.pi * position * ( 2.0 * j + 1.0 ) ) / ( 2 * self.height ) )
                summationResult += term_1 * term_2 * term_3
            # barP0 = self.alpha * self.Q / ( self.M + self.alpha * self.alpha * self.Q ) * \
            #     ( self.tao_0 + 0.5 * self.rho * self.g * self.height ) - self.rho_f * self.g * \
            #     ( 0.5 * self.height );
            barP0 = self._calculate_p_0(self.height)
            displacementValue = 8.0 * self.alpha * self.height * barP0 * summationResult / \
                ( math.pi * math.pi * self.M )
            displacementValue -= ( self.g / self.M ) * ( self.rho - self.alpha * self.rho_f ) * \
                ( self.height * yPosition - 0.5 * yPosition**2.0 )
            displacementValue -= ( self.tao_0 / self.M ) *yPosition
            return displacementValue;

    def __getPositionValuesAndSize( self, ny ):
        if type( ny ) == int:
            positionValues = self.getPositionValues( ny );
            size = ny
        elif type( ny ) == np.ndarray:
            positionValues = ny
            size = ny.size
        elif type( ny ) == list:
            positionValues = ny
            size = len( ny )
        return positionValues, size

    def getPressureValuesConstTime( self, time, numberOfSummationTerms = 200, ny = 200 ):
        '''If ny is an interger, then a list of equally spaced vertical position will be created. However, ny can be
            a list or an numpy array with specified vertical positions.'''
        positionValues, size = self.__getPositionValuesAndSize( ny )
        pressureValues = np.zeros(ny);
        for i in range( 0, size ):
            pressureValues[i] = self.getPressureValue( positionValues[ i ], time, numberOfSummationTerms )
            # pressureValue = self.getPressureValue( positionValues[ i ], time, numberOfSummationTerms );
            # pressureValues.append( pressureValue );
        # return np.array(pressureValues)
        return pressureValues

    def getDisplacementValuesConstTime( self, time, numberOfSummationTerms = 200, ny = 200 ):
        '''If ny is an interger, then a list of equally spaced vertical position will be created. However, ny can be
            a list or an numpy array with specified vertical positions.'''
        positionValues, size = self.__getPositionValuesAndSize( ny )
        displacementValues = [ ];
        for i in range( 0, size ):
            displacementValue = self.getDisplacementValue( positionValues[ i ], time, numberOfSummationTerms );
            displacementValues.append( displacementValue );
        return displacementValues;

    def getTimeValues( self, totalTimeInterval, timePoints = 200 ):
        timeValues = np.linspace(0, totalTimeInterval, timePoints )
        return timeValues;

    def getPressureValuesConstPosition( self, position, totalTimeInterval, numberOfSummationTerms = 200, timePoints = 200 ):
        timeValues = self.getTimeValues( totalTimeInterval, timePoints );
        pressureValues = [ ];
        for i in range( 0, len( timeValues ) ):
            pressureValue = self.getPressureValue( position, timeValues[ i ], numberOfSummationTerms );
            pressureValues.append( pressureValue );
        return pressureValues;

    def getDisplacementValuesConstPosition( self, position, totalTimeInterval, numberOfSummationTerms = 200, timePoints = 200 ):
        timeValues = self.getTimeValues( totalTimeInterval, timePoints );
        displacementValues = [ ];
        for i in range( 0, len( timeValues ) ):
            displacementValue = self.getDisplacementValue( position, timeValues[ i ], numberOfSummationTerms );
            displacementValues.append( displacementValue );
        return displacementValues;

def pEquilibrium( stress, height, c_f, c_s, phi, ni, G, alpha ):
    Lambda = 2*G*ni/(1-2*ni)
    Psi = phi*c_f + ( alpha - phi )*c_s
    Beta = Lambda*(1-ni)/ni + alpha*alpha/Psi
    p = -alpha*stress/( Psi*Beta )
    return p

def uEquilibrium( stress, height, c_f, c_s, phi, ni, G, alpha ):
    Psi = phi*c_f + ( alpha - phi )*c_s
    lame = 2*G*ni/(1-2*ni)
    Beta = lame*(1-ni)/ni + alpha*alpha/Psi
    return stress*height/Beta


