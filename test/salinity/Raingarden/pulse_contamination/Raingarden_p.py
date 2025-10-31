from proteus import *
from proteus.default_p import *
from proteus.richards import Richards
from domain_rg import *



analyticalSolution = None
    
viscosity     = 8.9e-4  #kg/(m*s)
density       = 998.2   #kg/m^3
gravity       = 9.8     #m/s^2
beta          = 0.0 #density*gravity*4.524e-10
m_per_s_by_m_per_d = 1.1574074e-5
lengthScale   = 1.0     #m
timeScale     = 1.0     #d #1.0/sqrt(g*lengthScale)
#make non-dimensional

permeability1  = (8.856*m_per_s_by_m_per_d)*viscosity/(gravity*density)  #m^2
#permeability1  = (0.00504)*viscosity/(gravity*density)  #m^2
thetaS1        = 0.40   #-
thetaR1        = 0.03   #-
mvg_alpha1     = 0.33   #1/m
mvg_n1         = 3.594
mvg_m1         = 1.0 - 1.0/mvg_n1
dimensionless_conductivity1  = (timeScale*density*gravity*permeability1/(viscosity*lengthScale))

permeability2  = (19.944*m_per_s_by_m_per_d)*viscosity/(gravity*density)  #m^2
#permeability2  = (0.00111)*viscosity/(gravity*density)  #m^2
thetaS2        = 0.37   #-
thetaR2        = 0.10   #-
mvg_alpha2     = 0.32   #1/m
mvg_n2         = 2.146
mvg_m2         = 1.0 - 1.0/mvg_n2
dimensionless_conductivity2  = (timeScale*density*gravity*permeability2/(viscosity*lengthScale))

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
    


class MyCoefficients(Richards.Coefficients):
    def attachModels(self, modelList):
        self.model = modelList[self.modelIndex]
        self.densityModel = modelList[1]  # TADR model
        self.model.q['rho']    = self.densityModel.q['rho']      
        self.model.ebqe['rho'] = self.densityModel.ebqe['rho']


LevelModelType = Richards.LevelModel
coefficients = Richards.Coefficients(nd,
                                     KsTypes,
                                     nVGtypes,
                                     alphaVGtypes,
                                     thetaRtypes,
                                     thetaSRtypes,
                                     gravity=dimensionless_gravity,
                                     density=dimensionless_density,
                                     beta=0.0001,
                                     DENSITY_MODEL=1,
                                     diagonal_conductivity=True,
                                     STABILIZATION_TYPE=2,
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

galerkin = False
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



pondingPressure=0.2

def getDBC_2D_Richards_Shock(x,flag):
    if x[1] == G[1]:
        if (x[0] >= G[0]/3.0 and
            x[0] <= 2.0*G[0]/3.0):
            return lambda x,t: pondingPressure
    if x[1] == 0.0:
        return lambda x,t: 0.0
    if (x[0] == 0.0 or
        x[0] == G[0]):
        return lambda x,t: x[1]*dimensionless_gravity[1]*dimensionless_density

dirichletConditions = {0:getDBC_2D_Richards_Shock}

class ShockIC_2D_Richards:
    def uOfXT(self,x,t):
        bc=getDBC_2D_Richards_Shock(x,0)
        if bc != None:
            return bc(x,t)
        else:
            return x[1]*dimensionless_gravity[1]*dimensionless_density

initialConditions  = {0:ShockIC_2D_Richards()}

fluxBoundaryConditions = {0:'noFlow'}

def getFBC_2D_Richards_Shock(x,flag):
    if x[1] == G[1]:
        if (x[0] < G[0]/3.0 or
            x[0] > 2.0*G[0]/3.0):
            return lambda x,t: 0.0

advectiveFluxBoundaryConditions =  {0:getFBC_2D_Richards_Shock}

diffusiveFluxBoundaryConditions = {0:{0:getFBC_2D_Richards_Shock}}

T = 100.0/timeScale
