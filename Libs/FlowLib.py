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
			linearSystem.addValueToVector(bIndex, -k*A/viscosity*density*gravity)
			linearSystem.addValueToVector(fIndex,  k*A/viscosity*density*gravity)

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
	# e = grid.getElements()[0]
	# r = grid.getRegions()[e.getParentRegionIndex()]
	# vIndex = e.getVertices()[0].getIndex()
	# linearSystem.addValueToMatrix(vIndex + pShift*n, vIndex + (1-pShift)*n, -biotOnRegions.getValue(r)/timeStep)

	# e = grid.getElements()[-1]
	# r = grid.getRegions()[e.getParentRegionIndex()]
	# vIndex = e.getVertices()[1].getIndex()
	# linearSystem.addValueToMatrix(vIndex + pShift*n, vIndex + (1-pShift)*n, biotOnRegions.getValue(r)/timeStep)

