from proteus import *
from proteus.default_p import *
from proteus.mphase_co2 import mphase_co2 

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

regions=[[0.1,0.1], [0.9, 1,1]]

regionFlags=[0,1]


domain = Domain.PlanarStraightLineGraphDomain(vertices= vertices,
                                              vertexFlags=vertexFlags,
                                              segments=segments,
                                              segmentFlags=segmentFlags,
                                              regions = regions,
                                              regionFlags = regionFlags,)
#dplt.plot_pslg_domain(polygon)

domain.writePoly('rg2d')
#if not regularGrid:
#    domain.writePoly('rg2d')
    #domain = Domain.PlanarStraightLineGraphDomain('rg2d')


analyticalSolution = None
    
viscosity     = 8.9e-4  #kg/(m*s)
density_water = 998.2   #kg/m^3
rho_air      = 1.225   #kg/m^3
gravity       = 9.8     #m/s^2
beta_water          = 0.0001 #density*gravity*4.524e-10
beta_air            = 0.0001 #density*gravity*4.524e-10
m_per_s_by_m_per_d = 1.1574074e-5
lengthScale   = 1.0     #m
timeScale     = 1.0     #d #1.0/sqrt(g*lengthScale)
#make non-dimensional

permeability1  = (3.438*m_per_s_by_m_per_d)*viscosity/(gravity*density_water)  #m^2
#permeability1  = (0.00504)*viscosity/(gravity*density)  #m^2
thetaS1        = 0.41   #-
thetaR1        = 0.02   #-
mvg_alpha1     = 1   #1/m
mvg_n1         = 1.6
mvg_m1         = 1.0 - 1.0/mvg_n1
dimensionless_conductivity1  = (timeScale*density_water*gravity*permeability1/(viscosity*lengthScale))

permeability2  = (1.60*m_per_s_by_m_per_d)*viscosity/(gravity*density_water)  #m^2
#permeability2  = (0.00111)*viscosity/(gravity*density)  #m^2
thetaS2        = 0.41   #-
thetaR2        = 0.041   #-
mvg_alpha2     = 1   #1/m
mvg_n2         = 1.378
mvg_m2         = 1.0 - 1.0/mvg_n2
dimensionless_conductivity2  = (timeScale*density_water*gravity*permeability2/(viscosity*lengthScale))

#print 'Ks',dimensionless_conductivity
dimensionless_density  = 1.0
dimensionless_gravity  = numpy.array([0.0,
                                        -1.0,
                                        0.0])
#dimensionless_alpha    = mvg_alpha*lengthScale
nMediaTypes  = 2
alphaVGtypes = numpy.zeros((nMediaTypes+1,),'d')
nVGtypes     = numpy.zeros((nMediaTypes+1,),'d')
thetaStypes  = numpy.zeros((nMediaTypes+1,),'d')
thetaRtypes  = numpy.zeros((nMediaTypes+1,),'d')
thetaSRtypes = numpy.zeros((nMediaTypes+1,),'d')
KsTypes      = numpy.zeros((nMediaTypes+1,2),'d')

for i in range(nMediaTypes+1):
    if i==0: 
        alphaVGtypes[i] = mvg_alpha1
        nVGtypes[i]     = mvg_n1
        thetaStypes[i]  = thetaS1
        thetaRtypes[i]  = thetaR1
        thetaSRtypes[i] = thetaStypes[i] - thetaRtypes[i]
        KsTypes[i,:]    = [dimensionless_conductivity1,dimensionless_conductivity1]#m/d?
    if i==1:
        alphaVGtypes[i] = mvg_alpha2
        nVGtypes[i]     = mvg_n2
        thetaStypes[i]  = thetaS2
        thetaRtypes[i]  = thetaR2
        thetaSRtypes[i] = thetaStypes[i] - thetaRtypes[i]
        KsTypes[i,:]    = [dimensionless_conductivity2,dimensionless_conductivity2]#m/d?    
    

LevelModelType = mphase_co2.LevelModel
coefficients = mphase_co2.Coefficients(nd,
                                     KsTypes,
                                     nVGtypes,
                                     alphaVGtypes,
                                     thetaRtypes,
                                     thetaSRtypes,
                                     gravity=dimensionless_gravity,
                                     density_water=dimensionless_density,
                                     density_air=rho_air/density_water,
                                     beta_water= beta_water, #=0.0001,
                                     beta_air = beta_air, #=0.0001,
                                     diagonal_conductivity=True,
                                     STABILIZATION_TYPE=0,
                                     ENTROPY_TYPE=1,
                                     LUMPED_MASS_MATRIX= False ,
                                     FCT=False,#True,
                                     num_fct_iter=0,
                                     # FOR ENTROPY VISCOSITY
                                     cE=1.0,
                                     uL=0.0,
                                     uR=1.0,
                                     # FOR ARTIFICIAL COMPRESSION
                                     cK=1.0,
                                     # OUTPUT quantDOFs
                                     outputQuantDOFs=False)
galerkin=False

#coefficients = ConservativeHeadRichardsMualemVanGenuchten(hydraulicConductivity=dimensionless_conductivity,
#                                                          gravity=dimensionless_gravity,
#                                                          density=dimensionless_density,
#                                                          thetaS=thetaS,
#                                                          thetaR=thetaR,
#                                                          alpha= dimensionless_alpha,
#                                                          n = mvg_n,
#                                                          m = mvg_m,
#                                                          beta = beta)


#coefficients = ConservativeHeadRichardsMualemVanGenuchten(hydraulicConductivity=dimensionless_conductivity,
#                                                          gravity=dimensionless_gravity,
#                                                          density=dimensionless_density,
#                                                          thetaS=thetaStypes,
#                                                          thetaR=thetaRtypes,
#                                                          alpha= alphaVGtypes,
#                                                          n = mvg_n,
#                                                          m = mvg_m,
#                                                          beta = beta)#



pondingPressure=0.5

def getDBC_2D_Richards_Shock_water(x,flag):
    if x[1] == G[1]:
        if (x[0] >= G[0]/3.0 and
            x[0] <= 2.0*G[0]/3.0):
            return lambda x,t: pondingPressure
    if x[1] == 0.0:
        return lambda x,t: 0.0
    if (x[0] == 0.0 or
        x[0] == G[0]):
        return lambda x,t: x[1]*dimensionless_gravity[1]*dimensionless_density

def getDBC_2D_Richards_Shock_air(x,flag):
    if x[1] == G[1]:
        if (x[0] >= G[0]/3.0 and
            x[0] <= 2.0*G[0]/3.0):
            return lambda x,t: 0.001
    if x[1] == 0.0:
        return lambda x,t: -0.001
    if (x[0] == 0.0 or
        x[0] == G[0]):
        return lambda x,t: -0.001 #x[1]*dimensionless_gravity[1]*dimensionless_density

dirichletConditions = {0:getDBC_2D_Richards_Shock_water, 1:getDBC_2D_Richards_Shock_air}

class ShockIC_2D_Richards_water:
    def uOfXT(self,x,t):
        bc=getDBC_2D_Richards_Shock_water(x,0)
        if bc != None:
            return bc(x,t)
        else:
            return x[1]*dimensionless_gravity[1]*dimensionless_density

class ShockIC_2D_Richards_air:
    def uOfXT(self,x,t):
        bc=getDBC_2D_Richards_Shock_air(x,0)
        if bc != None:
            return bc(x,t)
        else:
            return -0.0001
            #return x[1]*dimensionless_gravity[1]*dimensionless_density


initialConditions  = {0:ShockIC_2D_Richards_water(), 1:ShockIC_2D_Richards_air()}

fluxBoundaryConditions = {0:'noFlow', 1:'noFlow'}

def getFBC_2D_Richards_Shock(x,flag):
    if x[1] == G[1]:
        if (x[0] < G[0]/3.0 or
            x[0] > 2.0*G[0]/3.0):
            return lambda x,t: 0.0

advectiveFluxBoundaryConditions =  {0:getFBC_2D_Richards_Shock, 1:getFBC_2D_Richards_Shock}

diffusiveFluxBoundaryConditions = {0:{0:getFBC_2D_Richards_Shock}, 1:{1:getFBC_2D_Richards_Shock}}

T = 10.0/timeScale


