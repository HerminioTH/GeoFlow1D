def AssemblyDarcyVelocitiesToMatrix(linearSystem, grid, viscosity, permeability, pShift=0):
	for region in grid.getRegions():
		k = permeability.getValue(region)
		for elem in region.getElements():
			dx = elem.getLength()
			face = elem.getFace()
			A = face.getArea()
			backVertex = face.getBackwardVertex()
			forVertex = face.getForwardVertex()
			bIndex = backVertex.getIndex() + pShift*grid.getNumberOfVertices()
			fIndex = forVertex.getIndex() + pShift*grid.getNumberOfVertices()
			diffusiveOperator = [k*A/viscosity, -k*A/viscosity]
			for i,v in enumerate(elem.getVertices()):
				flux = diffusiveOperator[i]
				vIndex = v.getIndex() + pShift*grid.getNumberOfVertices()
				linearSystem.addValueToMatrix(bIndex, vIndex, +flux/dx)
				linearSystem.addValueToMatrix(fIndex, vIndex, -flux/dx)



def AssemblyDarcyVelocitiesToVector(linearSystem, grid, viscosity, permeability, density, gravity, pShift=0):
	n = grid.getNumberOfVertices()
	for region in grid.getRegions():
		k = permeability.getValue(region)
		for elem in region.getElements():
			dx = elem.getLength()
			face = elem.getFace()
			A = face.getArea()
			value = k*A*density*gravity/viscosity
			bIndex = face.getBackwardVertex().getIndex() + pShift*n
			fIndex = face.getForwardVertex().getIndex() + pShift*n
			linearSystem.addValueToVector(bIndex, -value)
			linearSystem.addValueToVector(fIndex,  value)


def AssemblyPisTermsToMatrix(linearSystem, grid, biotOnRegions, modulusOnRegions, timeStep, pShift=0):
	for region in grid.getRegions():
		M = modulusOnRegions.getValue(region)
		alpha = biotOnRegions.getValue(region)
		for elem in region.getElements():
			dx = elem.getLength()
			face = elem.getFace()
			A = face.getArea()
			backVertex = face.getBackwardVertex()
			forVertex = face.getForwardVertex()
			bIndex = backVertex.getIndex() + pShift*grid.getNumberOfVertices()
			fIndex = forVertex.getIndex() + pShift*grid.getNumberOfVertices()
			value = alpha*alpha*dx*A/(8*M*timeStep)
			diffusiveOperator = [value, -value]
			for i,v in enumerate(elem.getVertices()):
				flux = diffusiveOperator[i]
				vIndex = v.getIndex() + pShift*grid.getNumberOfVertices()
				linearSystem.addValueToMatrix(bIndex, vIndex, +flux/dx)
				linearSystem.addValueToMatrix(fIndex, vIndex, -flux/dx)


def AssemblyPisTermsToVector(linearSystem, biotOnRegions, modulusOnRegions, timeStep, p_old, pShift=0):
	n = grid.getNumberOfVertices()
	for region in grid.getRegions():
		M = modulusOnRegions.getValue(region)
		alpha = biotOnRegions.getValue(region)
		for e in region.getElements():
			dx = elem.getLength()
			face = elem.getFace()
			A = face.getArea()
			bVertex = face.getBackwardVertex()
			fVertex = face.getForwardVertex()
			pb = p_old.getValue(bVertex)
			pf = p_old.getValue(fVertex)
			value = -alpha*alpha*dx*A/(8*M*timeStep)
			linearSystem.addValueToVector(bVertex.getIndex() + pShift*n, value*(pf - pb))
			linearSystem.addValueToVector(fVertex.getIndex() + pShift*n, -value*(pf - pb))




def AssemblyFluidFlowAccumulationToMatrix(linearSystem, grid, timeStep, phiOnRegions, csOnRegions, cf, pShift=0):
	n = grid.getNumberOfVertices()
	for region in grid.getRegions():
		phi = phiOnRegions.getValue(region)
		cs = csOnRegions.getValue(region)
		for element in region.getElements():
			bIndex = element.getVertices()[0].getIndex()
			fIndex = element.getVertices()[1].getIndex()
			value = phi*(cs + cf)*element.getSubVolume()/timeStep
			linearSystem.addValueToMatrix(bIndex + pShift*n, bIndex + pShift*n, value)
			linearSystem.addValueToMatrix(fIndex + pShift*n, fIndex + pShift*n, value)

def AssemblyFluidFlowAccumulationToVector(linearSystem, grid, timeStep, phiOnRegions, csOnRegions, cf, p_old, pShift=0):
	n = grid.getNumberOfVertices()
	for region in grid.getRegions():
		phi = phiOnRegions.getValue(region)
		cs = csOnRegions.getValue(region)
		for element in region.getElements():
			bVertex = element.getVertices()[0]
			fVertex = element.getVertices()[1]
			value = phi*(cs + cf)*element.getSubVolume()/timeStep
			linearSystem.addValueToVector(bVertex.getIndex() + pShift*n, value*p_old.getValue(bVertex))
			linearSystem.addValueToVector(fVertex.getIndex() + pShift*n, value*p_old.getValue(fVertex))

def AssemblyBiotAccumulationToMatrix(linearSystem, grid, timeStep, biotOnRegions, phiOnRegions, csOnRegions, cf, pShift=0):
	n = grid.getNumberOfVertices()
	for region in grid.getRegions():
		phi = phiOnRegions.getValue(region)
		alpha = biotOnRegions.getValue(region)
		cs = csOnRegions.getValue(region)
		for element in region.getElements():
			bIndex = element.getVertices()[0].getIndex()
			fIndex = element.getVertices()[1].getIndex()
			value = (cf*phi + cs*(alpha - phi))*element.getSubVolume()/timeStep
			linearSystem.addValueToMatrix(bIndex + pShift*n, bIndex + pShift*n, value)
			linearSystem.addValueToMatrix(fIndex + pShift*n, fIndex + pShift*n, value)

def AssemblyBiotAccumulationToVector(linearSystem, grid, timeStep, biotOnRegions, phiOnRegions, csOnRegions, cf, p_old, pShift=0):
	n = grid.getNumberOfVertices()
	for region in grid.getRegions():
		phi = phiOnRegions.getValue(region)
		alpha = biotOnRegions.getValue(region)
		cs = csOnRegions.getValue(region)
		for element in region.getElements():
			bVertex = element.getVertices()[0]
			fVertex = element.getVertices()[1]
			value = (cf*phi + cs*(alpha - phi))*element.getSubVolume()/timeStep
			linearSystem.addValueToVector(bVertex.getIndex() + pShift*n, value*p_old.getValue(bVertex))
			linearSystem.addValueToVector(fVertex.getIndex() + pShift*n, value*p_old.getValue(fVertex))


def AssemblyFixedStressAccumulationToMatrix(linearSystem, grid, timeStep, biotOnRegions, deltaOnVertices, bulkModulusOnRegions, pShift=0):
	n = grid.getNumberOfVertices()
	for region in grid.getRegions():
		alpha = biotOnRegions.getValue(region)
		modulus = bulkModulusOnRegions.getValue(region)
		for element in region.getElements():
			bVertex = element.getVertices()[0]
			fVertex = element.getVertices()[1]
			bIndex = bVertex.getIndex()
			fIndex = fVertex.getIndex()
			value = (alpha*alpha/modulus)*element.getSubVolume()/timeStep
			linearSystem.addValueToMatrix(bIndex + pShift*n, bIndex + pShift*n, value/deltaOnVertices.getValue(bVertex))
			linearSystem.addValueToMatrix(fIndex + pShift*n, fIndex + pShift*n, value/deltaOnVertices.getValue(fVertex))

def AssemblyFixedStressAccumulationToVector(linearSystem, grid, timeStep, biotOnRegions, deltaOnVertices, bulkModulusOnRegions, p_old, pShift=0):
	n = grid.getNumberOfVertices()
	for region in grid.getRegions():
		alpha = biotOnRegions.getValue(region)
		modulus = bulkModulusOnRegions.getValue(region)
		for element in region.getElements():
			bVertex = element.getVertices()[0]
			fVertex = element.getVertices()[1]
			bIndex = bVertex.getIndex()
			fIndex = fVertex.getIndex()
			value = (alpha*alpha/modulus)*element.getSubVolume()/timeStep
			linearSystem.addValueToVector(bVertex.getIndex() + pShift*n, value*p_old.getValue(bVertex)/deltaOnVertices.getValue(bVertex))
			linearSystem.addValueToVector(fVertex.getIndex() + pShift*n, value*p_old.getValue(fVertex)/deltaOnVertices.getValue(fVertex))








def AssemblyFixedStressAccumulationToMatrix_2(linearSystem, grid, timeStep, biotOnRegions, deltaOnRegions, bulkModulusOnRegions, pShift=0):
	n = grid.getNumberOfVertices()
	for region in grid.getRegions():
		alpha = biotOnRegions.getValue(region)
		modulus = bulkModulusOnRegions.getValue(region)
		delta = deltaOnRegions.getValue(region)
		for element in region.getElements():
			bVertex = element.getVertices()[0]
			fVertex = element.getVertices()[1]
			bIndex = bVertex.getIndex()
			fIndex = fVertex.getIndex()
			value = (alpha*alpha/modulus)*element.getSubVolume()/timeStep/delta
			linearSystem.addValueToMatrix(bIndex + pShift*n, bIndex + pShift*n, value)
			linearSystem.addValueToMatrix(fIndex + pShift*n, fIndex + pShift*n, value)

def AssemblyFixedStressAccumulationToVector_2(linearSystem, grid, timeStep, biotOnRegions, deltaOnRegions, bulkModulusOnRegions, p_old, pShift=0):
	n = grid.getNumberOfVertices()
	for region in grid.getRegions():
		alpha = biotOnRegions.getValue(region)
		modulus = bulkModulusOnRegions.getValue(region)
		delta = deltaOnRegions.getValue(region)
		for element in region.getElements():
			bVertex = element.getVertices()[0]
			fVertex = element.getVertices()[1]
			bIndex = bVertex.getIndex()
			fIndex = fVertex.getIndex()
			value = (alpha*alpha/modulus)*element.getSubVolume()/timeStep/delta
			linearSystem.addValueToVector(bVertex.getIndex() + pShift*n, value*p_old.getValue(bVertex))
			linearSystem.addValueToVector(fVertex.getIndex() + pShift*n, value*p_old.getValue(fVertex))


def AssemblyVolumetricStrainToMatrix(linearSystem, grid, timeStep, biotOnRegions, pShift=0):
	n = grid.getNumberOfVertices()
	for region in grid.getRegions():
		value = biotOnRegions.getValue(region)/(2*timeStep)
		for e in region.getElements():
			bVertex = e.getVertices()[0]
			fVertex = e.getVertices()[1]
			linearSystem.addValueToMatrix(bVertex.getIndex() + pShift*n, bVertex.getIndex() + (1-pShift)*n, value)
			linearSystem.addValueToMatrix(bVertex.getIndex() + pShift*n, fVertex.getIndex() + (1-pShift)*n, value)
			linearSystem.addValueToMatrix(fVertex.getIndex() + pShift*n, bVertex.getIndex() + (1-pShift)*n, -value)
			linearSystem.addValueToMatrix(fVertex.getIndex() + pShift*n, fVertex.getIndex() + (1-pShift)*n, -value)
	e = grid.getElements()[0]
	r = grid.getRegions()[e.getParentRegionIndex()]
	vIndex = e.getVertices()[0].getIndex()
	linearSystem.addValueToMatrix(vIndex + pShift*n, vIndex + (1-pShift)*n, -biotOnRegions.getValue(r)/timeStep)

	e = grid.getElements()[-1]
	r = grid.getRegions()[e.getParentRegionIndex()]
	vIndex = e.getVertices()[1].getIndex()
	linearSystem.addValueToMatrix(vIndex + pShift*n, vIndex + (1-pShift)*n, biotOnRegions.getValue(r)/timeStep)

def AssemblyVolumetricStrainToVector(linearSystem, grid, timeStep, biotOnRegions, u_old, pShift=0):
	n = grid.getNumberOfVertices()
	for region in grid.getRegions():
		value = biotOnRegions.getValue(region)/(2*timeStep)
		for e in region.getElements():
			bVertex = e.getVertices()[0]
			fVertex = e.getVertices()[1]
			ub = u_old.getValue(bVertex)
			uf = u_old.getValue(fVertex)
			linearSystem.addValueToVector(bVertex.getIndex() + pShift*n, value*(ub + uf))
			linearSystem.addValueToVector(fVertex.getIndex() + pShift*n, -value*(ub + uf))
	e = grid.getElements()[0]
	r = grid.getRegions()[e.getParentRegionIndex()]
	bVertex = e.getVertices()[0]
	linearSystem.addValueToVector(bVertex.getIndex() + pShift*n, -biotOnRegions.getValue(r)*u_old.getValue(bVertex)/(1*timeStep))

	e = grid.getElements()[-1]
	r = grid.getRegions()[e.getParentRegionIndex()]
	fVertex = e.getVertices()[1]
	linearSystem.addValueToVector(fVertex.getIndex() + pShift*n, biotOnRegions.getValue(r)*u_old.getValue(fVertex)/(1*timeStep))



def AssemblyVolumetricStrainToVector2(linearSystem, grid, timeStep, biotOnRegions, u_new, u_old, pShift=0):
	n = grid.getNumberOfVertices()
	for region in grid.getRegions():
		value = biotOnRegions.getValue(region)/(2*timeStep)
		for e in region.getElements():
			bVertex = e.getVertices()[0]
			fVertex = e.getVertices()[1]
			uob = u_old.getValue(bVertex)
			uof = u_old.getValue(fVertex)
			unb = u_new.getValue(bVertex)
			unf = u_new.getValue(fVertex)
			linearSystem.addValueToVector(bVertex.getIndex() + pShift*n, value*(uob + uof - unb - unf))
			linearSystem.addValueToVector(fVertex.getIndex() + pShift*n, -value*(uob + uof - unb - unf))
	e = grid.getElements()[0]
	r = grid.getRegions()[e.getParentRegionIndex()]
	bVertex = e.getVertices()[0]
	linearSystem.addValueToVector(bVertex.getIndex() + pShift*n, -biotOnRegions.getValue(r)*(u_old.getValue(bVertex) - u_new.getValue(bVertex))/timeStep)

	e = grid.getElements()[-1]
	r = grid.getRegions()[e.getParentRegionIndex()]
	fVertex = e.getVertices()[1]
	linearSystem.addValueToVector(fVertex.getIndex() + pShift*n, biotOnRegions.getValue(r)*(u_old.getValue(fVertex) - u_new.getValue(fVertex))/timeStep)






# ---------------------------- LOOP BY ELEMENTS ----------------------------------


def AssemblyBiotAccumulationToMatrix_e(linearSystem, grid, timeStep, biotOnElements, phiOnElements, csOnElements, cf, pShift=0):
	n = grid.getNumberOfVertices()
	for element in grid.getElements():
		phi = phiOnElements.getValue(element)
		alpha = biotOnElements.getValue(element)
		cs = csOnElements.getValue(element)
		bIndex = element.getVertices()[0].getIndex()
		fIndex = element.getVertices()[1].getIndex()
		value = (cf*phi + cs*(alpha - phi))*element.getSubVolume()/timeStep
		linearSystem.addValueToMatrix(bIndex + pShift*n, bIndex + pShift*n, value)
		linearSystem.addValueToMatrix(fIndex + pShift*n, fIndex + pShift*n, value)

def AssemblyBiotAccumulationToVector_e(linearSystem, grid, timeStep, biotOnElements, phiOnElements, csOnElements, cf, p_old, pShift=0):
	n = grid.getNumberOfVertices()
	for element in grid.getElements():
		phi = phiOnElements.getValue(element)
		alpha = biotOnElements.getValue(element)
		cs = csOnElements.getValue(element)
		bVertex = element.getVertices()[0]
		fVertex = element.getVertices()[1]
		value = (cf*phi + cs*(alpha - phi))*element.getSubVolume()/timeStep
		linearSystem.addValueToVector(bVertex.getIndex() + pShift*n, value*p_old.getValue(bVertex))
		linearSystem.addValueToVector(fVertex.getIndex() + pShift*n, value*p_old.getValue(fVertex))

def AssemblyFixedStressAccumulationToMatrix_e(linearSystem, grid, timeStep, biotOnElements, deltaOnElements, bulkModulusOnElements, pShift=0):
	n = grid.getNumberOfVertices()
	for element in grid.getElements():
		alpha = biotOnElements.getValue(element)
		modulus = bulkModulusOnElements.getValue(element)
		delta = deltaOnElements.getValue(element)
		bVertex = element.getVertices()[0]
		fVertex = element.getVertices()[1]
		bIndex = bVertex.getIndex()
		fIndex = fVertex.getIndex()
		value = (alpha*alpha/modulus)*element.getSubVolume()/timeStep/delta
		linearSystem.addValueToMatrix(bIndex + pShift*n, bIndex + pShift*n, value)
		linearSystem.addValueToMatrix(fIndex + pShift*n, fIndex + pShift*n, value)

def AssemblyFixedStressAccumulationToVector_e(linearSystem, grid, timeStep, biotOnElements, deltaOnElements, bulkModulusOnElements, p_old, pShift=0):
	n = grid.getNumberOfVertices()
	for element in grid.getElements():
		alpha = biotOnElements.getValue(element)
		modulus = bulkModulusOnElements.getValue(element)
		delta = deltaOnElements.getValue(element)
		bVertex = element.getVertices()[0]
		fVertex = element.getVertices()[1]
		bIndex = bVertex.getIndex()
		fIndex = fVertex.getIndex()
		value = (alpha*alpha/modulus)*element.getSubVolume()/timeStep/delta
		linearSystem.addValueToVector(bVertex.getIndex() + pShift*n, value*p_old.getValue(bVertex))
		linearSystem.addValueToVector(fVertex.getIndex() + pShift*n, value*p_old.getValue(fVertex))

def AssemblyDarcyVelocitiesToMatrix_e(linearSystem, grid, viscosity, permeabilityOnElements, pShift=0):
	for element in grid.getElements():
		k = permeabilityOnElements.getValue(element)
		dx = element.getLength()
		face = element.getFace()
		A = face.getArea()
		backVertex = face.getBackwardVertex()
		forVertex = face.getForwardVertex()
		bIndex = backVertex.getIndex() + pShift*grid.getNumberOfVertices()
		fIndex = forVertex.getIndex() + pShift*grid.getNumberOfVertices()
		diffusiveOperator = [k*A/viscosity, -k*A/viscosity]
		for i,v in enumerate(element.getVertices()):
			flux = diffusiveOperator[i]
			vIndex = v.getIndex() + pShift*grid.getNumberOfVertices()
			linearSystem.addValueToMatrix(bIndex, vIndex, +flux/dx)
			linearSystem.addValueToMatrix(fIndex, vIndex, -flux/dx)

def AssemblyDarcyVelocitiesToVector_e(linearSystem, grid, viscosity, permeabilityOnElements, density, gravity, pShift=0):
	n = grid.getNumberOfVertices()
	for element in grid.getElements():
		k = permeabilityOnElements.getValue(element)
		dx = element.getLength()
		face = element.getFace()
		A = face.getArea()
		value = k*A*density*gravity/viscosity
		bIndex = face.getBackwardVertex().getIndex() + pShift*n
		fIndex = face.getForwardVertex().getIndex() + pShift*n
		linearSystem.addValueToVector(bIndex, -value)
		linearSystem.addValueToVector(fIndex,  value)

def AssemblyVolumetricStrainToVector_e(linearSystem, grid, timeStep, biotOnElements, u_old, pShift=0):
	n = grid.getNumberOfVertices()
	for element in grid.getElements():
		value = biotOnElements.getValue(element)/(2*timeStep)
		bVertex = element.getVertices()[0]
		fVertex = element.getVertices()[1]
		ub = u_old.getValue(bVertex)
		uf = u_old.getValue(fVertex)
		linearSystem.addValueToVector(bVertex.getIndex() + pShift*n, value*(ub + uf))
		linearSystem.addValueToVector(fVertex.getIndex() + pShift*n, -value*(ub + uf))
	e = grid.getElements()[0]
	bVertex = e.getVertices()[0]
	linearSystem.addValueToVector(bVertex.getIndex() + pShift*n, -biotOnElements.getValue(e)*u_old.getValue(bVertex)/(1*timeStep))

	e = grid.getElements()[-1]
	fVertex = e.getVertices()[1]
	linearSystem.addValueToVector(fVertex.getIndex() + pShift*n, biotOnElements.getValue(e)*u_old.getValue(fVertex)/(1*timeStep))
