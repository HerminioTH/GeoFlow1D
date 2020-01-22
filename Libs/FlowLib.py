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
			diffusiveOperator = [k*A/viscosity, -k*A/viscosity]
			linearSystem.addValueToVector(bIndex, -diffusiveOperator[0]*density*gravity)
			linearSystem.addValueToVector(fIndex, -diffusiveOperator[1]*density*gravity)



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