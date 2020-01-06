import unittest

def AssemblyStiffnessMatrix(linearSystem, grid, props, uShift):
	for e in grid.getElements():
		dx = e.getLength()
        f = e.getFace()
        A = f.getArea()
        backVertex = f.getBackwardVertex()
        forVertex = f.getForwardVertex()
        bIndex = backVertex.getIndex() + uShift*grid.getNumberOfVertices()
        fIndex = forVertex.getIndex() + uShift*grid.getNumberOfVertices()
        value = props.lame.getValue( e )*(1-props.poisson)/props.poisson
        forceOperator = [-value/dx, value/dx]
        localIndex = 0
        for v in e.getVertices():
            flux = forceOperator[localIndex]
            vIndex = v.getIndex() + uShift*grid.getNumberOfVertices()
            linearSystem.addValueToMatrix( bIndex, vIndex, flux )
            linearSystem.addValueToMatrix( fIndex, vIndex, -flux )
            localIndex += 1

def AssemblyPorePressureToGeoMatrix(linearSystem, grid, props, uShift):
	for e in grid.getElements():
		f = e.getFace()
        A = f.getArea()
        backVertex = f.getBackwardVertex()
        forVertex = f.getForwardVertex()
        bIndex = backVertex.getIndex() + uShift*grid.getNumberOfVertices()
        fIndex = forVertex.getIndex() + uShift*grid.getNumberOfVertices()
        operator = [props.biot*A/2., props.biot*A/2.]
        localIndex = 0
        for v in e.getVertices():
            flux = operator[localIndex]
            vIndex = v.getIndex() + (1-uShift)*grid.getNumberOfVertices()
            ls.addValueToMatrix( bIndex, vIndex, flux )
            ls.addValueToMatrix( fIndex, vIndex, -flux )
            localIndex += 1
	ls.addValueToMatrix( 2*nv-1, nv-1, 2*biot*A )