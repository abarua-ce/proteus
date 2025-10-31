from proteus import *
from proteus.default_p import *
from proteus.richards import Richards

nd = 2

#L=(10.0,10.0,1.0)

G=(1.8, 1.2,1)
regularGrid=False

domain = Domain.PlanarStraightLineGraphDomain()
boundaries=['bottom','top','left','right',
            'Storage','Root']

boundaryTags=dict([(key,i+1) for (i,key) in enumerate(boundaries)])

vertices= [[0.0,0.0], #0
           [G[0] ,0.0], #1
           [G[0] ,0.7], #2
           [G[0], G[1]], #3
           [0.0, G[1]], #4 
           [0.0, 0.7], #5
          ]

vertexFlags=[boundaryTags['bottom'],
             boundaryTags['bottom'],
             boundaryTags['right'],
             boundaryTags['top'],
             boundaryTags['top'],
             boundaryTags['left'],
             ]


segments=[[0,1],
          [1,2],
          [2,5],
          [2,3],
          [3,4],
          [4,5],
          [5,0]]

segmentFlags=[boundaryTags['bottom'],
              boundaryTags['right'],
              boundaryTags['Storage'],
              boundaryTags['right'],
              boundaryTags['top'],
              boundaryTags['left'],
              boundaryTags['left']]

regions=[[0.1,0.1], [0.9, 1.1]]

regionFlags=[0,1]


domain = Domain.PlanarStraightLineGraphDomain(vertices= vertices,
                                              vertexFlags=vertexFlags,
                                              segments=segments,
                                              segmentFlags=segmentFlags,
                                              regions = regions,
                                              regionFlags = regionFlags,)
#dplt.plot_pslg_domain(polygon)


if not regularGrid:
    domain.writePoly('rg2d')
    #domain = Domain.PlanarStraightLineGraphDomain('rg2d')

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}
#femSpaces = {0:C0_AffineQuadraticOnSimplexWithNodalBasis}

#elementQuadrature = SimplexGaussQuadrature(nd,4)

#elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,4)

elementQuadrature = SimplexLobattoQuadrature(nd,1)
#
elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,1)

nnx=41
nny=41
nLevels = 1
triangleFlag = 0
triangleOptions="pAq30Dena%f" % (0.5*(L[0]/(nnx-1))**2,)
#he= 0.4

#triangleOptions="pa0.02"


subgridError = None
