def AssemblyStiffnessMatrix(linearSystem, grid, modulus, uShift):
    for region in grid.getRegions():
        value = modulus.getValue(region)
    	for e in region.getElements():
            dx = e.getLength()
            f = e.getFace()
            bIndex = f.getBackwardVertex().getIndex() + uShift*grid.getNumberOfVertices()
            fIndex = f.getForwardVertex().getIndex() + uShift*grid.getNumberOfVertices()
            forceOperator = [-value/dx, value/dx]
            localIndex = 0
            for v in e.getVertices():
                flux = forceOperator[localIndex]
                vIndex = v.getIndex() + uShift*grid.getNumberOfVertices()
                linearSystem.addValueToMatrix( bIndex, vIndex, flux )
                linearSystem.addValueToMatrix( fIndex, vIndex, -flux )
                localIndex += 1

def AssemblyGravityToVector(linearSystem, grid, density, gravity, uShift):
    n = grid.getNumberOfVertices()
    for region in grid.getRegions():
        rho = density.getValue(region)
        for elem in region.getElements():
            face = elem.getFace()
            bVertex = face.getBackwardVertex()
            fVertex = face.getForwardVertex()
            value = -rho*gravity*elem.getSubVolume()
            linearSystem.addValueToVector(bVertex.getIndex() + uShift*n, value)
            linearSystem.addValueToVector(fVertex.getIndex() + uShift*n, value)

def AssemblyPorePressureToMatrix(linearSystem, grid, biot, uShift):
    for region in grid.getRegions():
        alpha = biot.getValue(region)
        for e in region.getElements():
            f = e.getFace()
            backVertex = f.getBackwardVertex()
            forVertex = f.getForwardVertex()
            bIndex = backVertex.getIndex() + uShift*grid.getNumberOfVertices()
            fIndex = forVertex.getIndex() + uShift*grid.getNumberOfVertices()
            operator = [-alpha/2., -alpha/2.]
            for i,v in enumerate(e.getVertices()):
                flux = operator[i]
                col = v.getIndex() + (1-uShift)*grid.getNumberOfVertices()
                linearSystem.addValueToMatrix( bIndex, col, flux )
                linearSystem.addValueToMatrix( fIndex, col, -flux )


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



if __name__ == '__main__':
    from GridLib import *
    from FieldsLib import *
    from LinearSystem import *

    L_0 = 4.
    L_1 = 6.
    L = L_0 + L_1
    nVertices = 10
    nodesCoord, elemConn = createGridData( L, nVertices )

    # -------------- GRID DATA ----------------------------
    gridData = GridData()
    gridData.setElementConnectivity( elemConn )
    gridData.setNodeCoordinates( nodesCoord )
    centroidCoord = []
    for e in elemConn:
        x_0 = gridData.nodeCoordinates[e[0]]
        x_1 = gridData.nodeCoordinates[e[1]]
        centroidCoord.append((x_0 + x_1)/2.)
    region_1 = []
    region_2 = []
    namesOfRegions = ['bottom', 'top']
    for e, x in enumerate(centroidCoord):
        if x <= L_0:
            region_1.append(e)
        elif x > L_0:
            region_2.append(e)
    gridData.setElementsToRegion(region_1, 'lower_layer')
    gridData.setElementsToRegion(region_2, 'upper_layer')
    g = Grid_1D( gridData )

    for region in g.getRegions():
        print region.getName()
        for element in region.getElements():
            vec = [element.getIndex()]
            for v in element.getVertices():
                vec.append(v.getIndex())
            print vec
        print '\n'
    # -----------------------------------------------------

    # -------------- PROPERTIES ----------------------------
    M = ScalarField(g.getNumberOfRegions())
    M.setValue(g.getRegions()[0], 1000.)
    M.setValue(g.getRegions()[1], 2000.)
    # -----------------------------------------------------

    # -------------- LINEAR SYSTEM ------------------------
    ls = LinearSystem(g.getNumberOfVertices())
    AssemblyStiffnessMatrix(ls, g, M, 0)
    ls.applyDirichlet(0, 0)
    ls.applyNeumann(-1, -1000)
    print g.getNumberOfVertices()
    print ls.getMatrix()
    print ls.getVector()
    ls.solve()
    print ls.getSolution()
    # -----------------------------------------------------


