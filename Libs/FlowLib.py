def AssemblyMassDarcyVelocities(linearSystem, grid, viscosity, permeability, density=0, gravity=0, pShift=0):
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
				vIndex = v.getIndex()
				linearSystem.addValueToMatrix(bIndex, vIndex, +flux/dx)
				linearSystem.addValueToMatrix(fIndex, vIndex, -flux/dx)
			linearSystem.addValueToVector(bIndex, -diffusiveOperator[0]*density*gravity)
			linearSystem.addValueToVector(fIndex, -diffusiveOperator[1]*density*gravity)



# def Assembly



