import unittest

def AssemblyStiffnessMatrix(linearSystem, grid, lame, poisson, uShift):
    for region in grid.getRegions():
        value = lame.getValue(region)*(1 - poisson.getValue(region))/poisson.getValue(region)
    	for e in region.getElements():
            dx = e.getLength()
            f = e.getFace()
            A = f.getArea()
            backVertex = f.getBackwardVertex()
            forVertex = f.getForwardVertex()
            bIndex = backVertex.getIndex() + uShift*grid.getNumberOfVertices()
            fIndex = forVertex.getIndex() + uShift*grid.getNumberOfVertices()
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
    R1 = []
    R2 = []
    namesOfRegions = ['bottom', 'top']
    for e, x in enumerate(centroidCoord):
        if x <= L_0:
            R1.append(e)
        elif x > L_0:
            R2.append(e)
    elemOnRegion1 = gridData.elemConnectivity[R1[0]:R1[-1]+1]
    elemOnRegion2 = gridData.elemConnectivity[R2[0]:R2[-1]+1]
    print gridData.elemConnectivity
    gridData.setElementsToRegions([elemOnRegion1, elemOnRegion2], namesOfRegions)
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
    lame = ScalarField(g.getNumberOfRegions())
    poisson = ScalarField(g.getNumberOfRegions())
    valuesLame = [1, 1]
    valuesPoisson = [0.1, 0.1]
    i = 0
    for region in g.getRegions():
        lame.setValue(region, valuesLame[i])
        poisson.setValue(region, valuesPoisson[i])
        i += 1
    # -----------------------------------------------------

    # -------------- LINEAR SYSTEM ------------------------
    ls = LinearSystem(g.getNumberOfVertices())
    AssemblyStiffnessMatrix(ls, g, lame, poisson, 0)
    ls.applyDirichlet(0, 0)
    ls.applyNeumann(-1, 1000)
    print g.getNumberOfVertices()
    print ls.getMatrix()
    print ls.getVector()
    ls.solve()
    print ls.getSolution()
    # -----------------------------------------------------


