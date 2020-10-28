class FrozenClass(object):
    __isfrozen = False
    def __setattr__(self, key, value):
        if self.__isfrozen and not hasattr(self, key):
            raise TypeError( "%r is a frozen class" % self )
        object.__setattr__(self, key, value)

    def _freeze(self):
        self.__isfrozen = True

class PhysicalProperties(FrozenClass):
	def __init__(self):
		self.k = None			# Permeability
		self.phi = None 		# Porosity
		self.c_s = None 		# Solid compressibility
		self.M = None 			# Constrained modulus
		self.biot = None		# Biot's coefficient
		self.rho = None 		# rho_f*phi + rho_s*(1 - phi)
		self.rho_f = None 		# Fluid density
		self.c_f = None 		# Fluid compressibility
		self.mu = None 			# Fluid viscosity

		self._freeze()
