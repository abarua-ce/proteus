from proteus import *
from proteus.default_p import *
from proteus.mphase_co2 import mphase_co2

nd = 1

L=(1.0,1.0,1.0)

analyticalSolution = None

viscosity     = 8.9e-4  #kg/(m*s)
density_water = 998.2   #kg/m^3
density_air   = 1.2
gravity       = 9.8     #m/s^2
beta_water    = 0.0
beta_air      = 0.0#density*gravity*4.524e-10

m_per_s_by_m_per_d = 1.1574074e-5
permeability  = (0.00922*m_per_s_by_m_per_d)*viscosity/(gravity*density_water)  #m^2
thetaS        = 0.368   #-
thetaR        = 0.102   #-
mvg_alpha     = 0.0335    #1/m
mvg_n         = 1.694
mvg_m         = 1.0 - 1.0/mvg_n
lengthScale   = 1.0     #m
timeScale     = 1.0     #d #1.0/sqrt(g*lengthScale)
#make non-dimensional
dimensionless_conductivity  = (timeScale*density_water*gravity*permeability/(viscosity*lengthScale))/m_per_s_by_m_per_d
dimensionless_density_water  = 1.0
dimensionless_density_air  = density_air/density_water

dimensionless_gravity  = numpy.array([-1.0,
                                       0.0,
                                       0.0])
dimensionless_alpha    = mvg_alpha*lengthScale
satRichards = False
optRichards = True
nMediaTypes  = 1
alphaVGtypes = numpy.zeros((nMediaTypes+1,),'d')
nVGtypes     = numpy.zeros((nMediaTypes+1,),'d')
thetaStypes  = numpy.zeros((nMediaTypes+1,),'d')
thetaRtypes  = numpy.zeros((nMediaTypes+1,),'d')
thetaSRtypes = numpy.zeros((nMediaTypes+1,),'d')
KsTypes      = numpy.zeros((nMediaTypes+1,1),'d')
BC_entry_head= numpy.zeros((nMediaTypes+1,),'d')
BC_lambda    = numpy.zeros((nMediaTypes+1,),'d')

for i in range(nMediaTypes+1):
    alphaVGtypes[i]     = mvg_alpha
    nVGtypes[i]         = mvg_n
    thetaStypes[i]      = thetaS
    thetaRtypes[i]      = thetaR
    thetaSRtypes[i]     = thetaStypes[i] - thetaRtypes[i]
    KsTypes[i,:]        = [dimensionless_conductivity]#,dimensionless_conductivity,dimensionless_conductivity]#m/d?
    BC_entry_head[i]  = 0.0
    BC_lambda[i]      = 2.0
	
useSeepageFace = True
galerkin=False

# if galerkin:
#     stabilization_type=0
# else:
#     stabilization_type=1

LevelModelType = mphase_co2.LevelModel
coefficients = mphase_co2.Coefficients(nd,
                                     KsTypes,
                                     nVGtypes,
                                     alphaVGtypes,
                                     thetaRtypes,
                                     thetaSRtypes,
                                     gravity=dimensionless_gravity,
                                     density_water =dimensionless_density_water,
                                     density_air = dimensionless_density_air,
                                     beta_water=beta_water,
                                     beta_air = beta_air,
                                     diagonal_conductivity=True,
                                     STABILIZATION_TYPE=2,
                                     PSK_TYPE =0,
                                     ENTROPY_TYPE=1,
                                     LUMPED_MASS_MATRIX=False,
                                     FCT=False,
                                     MONOLITHIC=False,
                                     num_fct_iter=1,
                                         # FOR ENTROPY VISCOSITY
                                     cE=1.0,
                                     uL=0.0,
                                     uR=1.0,
                                     # FOR ARTIFICIAL COMPRESSION
                                     cK=1.0,
                                     # OUTPUT quantDOFs
                                     outputQuantDOFs=False, 
                                     BC_entry_head= BC_entry_head,
                                     BC_lambda = BC_lambda)


#pondingPressure=-0.1#-0.1
#bottomPressure = -0.2#0.0
pondingPressure= -0.75 #0.1
bottomPressure = -10.0
#pondingSaturation = 0.9
#waterTableSaturation = 0.9
#initialSaturation = 0.01
#pondingPressure=-0.1
# if satRichards:
#     def getDBC_Richards_Shock(x,flag):
#         if x[0] == L[0]:
#             return lambda x,t: pondingSaturation
#         if x[0] == 0.0:
#             return lambda x,t: waterTableSaturation
#else:
def getDBC_Richards_Shock_water(x,flag):
    if x[0] == L[0]:
        return lambda x,t: pondingPressure
    if x[0] == 0.0:
        return lambda x,t: bottomPressure
   
def getDBC_Richards_Shock_air(x,flag):
    if x[0] == L[0]:
        return lambda x,t: pondingPressure
    if x[0] == 0.0:
        return lambda x,t: bottomPressure


dirichletConditions = {0:getDBC_Richards_Shock_water,
                       1:getDBC_Richards_Shock_air}

# if satRichards:
#     class ShockIC_Richards:
#         def uOfXT(self,x,t):
#             f = getDBC_Richards_Shock(x,0)
#             if f:
#                 return f(x,t)
#             return initialSaturation
# else:
class ShockIC_Richards_water:
    def uOfXT(self,x,t):
        f = getDBC_Richards_Shock_water(x,0)
        if f:
            return f(x,t)
        else:
            return -10.0
        #     # return bottomPressure + x[0]*dimensionless_gravity[0]*dimensionless_density
        # if x[0] < L[0]:#*0.5:
        #     return bottomPressure + x[0]*dimensionless_gravity[0]*dimensionless_density
        # else:
        #     return pondingPressure
class ShockIC_Richards_air:
    def uOfXT(self,x,t):
        f = getDBC_Richards_Shock_air(x,0)
        if f:
            return f(x,t)
        else:
            return -10.0


initialConditions  = {0:ShockIC_Richards_water(),
                      1:ShockIC_Richards_air()}

fluxBoundaryConditions = {0:'outFlow',
                          1: 'outFlow'}

def flux(x,flag):
    return None
#    if x[0] == L[0]:
#        return lambda x,t: 0.0
#    if x[0] == 0.0:
#        return lambda x,t: 0.0

advectiveFluxBoundaryConditions =  {0:flux, 1:flux}

diffusiveFluxBoundaryConditions = {0:{}, 1:{}}

T = 1/24/timeScale
#T = 0.35/timeScale
