from proteus import *
from proteus.default_p import *
from proteus.richards import Richards


import math
nd = 3

regularGrid=False

boundaries = ['outer', 'top', 'bottom', 'pump']
boundaryTags = {key: i + 1 for i, key in enumerate(boundaries)}


def generate_hollow_cylinder(center, inner_radius=0.1, outer_radius=30,
                              height=3.0, n_points=20):
    #cx, cy = center
    cx, cy, cz = center
    angle_step = 2 * math.pi / n_points

    outer_top = []
    outer_bottom = []
    inner_top = []
    inner_bottom = []

    for i in range(n_points):
        angle = i * angle_step
        dx = math.cos(angle)
        dy = math.sin(angle)
        outer_top.append([cx + outer_radius * dx, cy + outer_radius * dy, cz + height/2])
        outer_bottom.append([cx + outer_radius * dx, cy + outer_radius * dy, cz-height/2])
        inner_top.append([cx + inner_radius * dx, cy + inner_radius * dy, cz + height/2])
        inner_bottom.append([cx + inner_radius * dx, cy + inner_radius * dy, cz- height/2])
    
    facets = []
    for i in range(n_points):
        ni = (i + 1) % n_points
        facets.append([i, ni, n_points + ni, n_points + i])  # outer wall
        facets.append([2 * n_points + ni, 2 * n_points + i, 3 * n_points + i, 3 * n_points + ni])  # inner wall
        facets.append([i, ni, 2 * n_points + ni, 2 * n_points + i])  # top face
        facets.append([3 * n_points + i, 3 * n_points + ni, n_points + ni, n_points + i])  # bottom face

    #return vertices, facets

    vertices = outer_top + outer_bottom + inner_top + inner_bottom
    facets_proteus = [[f] for f in facets]
    return vertices, facets_proteus



center=[30.0, 30.0, 1.5]
pump_radius= 0.6
domain_radius= 30.0
height= 3.0
n_points= 10

vertices, facets = generate_hollow_cylinder(center, pump_radius, domain_radius, height,n_points)


vertexFlags = (
    [boundaryTags['top']] * n_points +   # outer_top
    [boundaryTags['bottom']] * n_points +   # outer_bottom
    [boundaryTags['pump']]  * n_points +   # inner_top
    [boundaryTags['pump']]  * n_points     # inner_bottom
    )

facetFlags = (
    [boundaryTags['outer']]  * n_points +  # outer wall
    [boundaryTags['pump']]   * n_points +  # inner wall
    [boundaryTags['top']]    * n_points +  # top ring face
    [boundaryTags['bottom']] * n_points    # bottom ring face
)

region_radius = (pump_radius + domain_radius) / 2.0  # (0.6 + 30) / 2 = 15.3
regions = [[center[0] + region_radius, center[1], center[2]]]

regionFlags= [1]

holes = [[30.0, 30.0, 1.5]]

domain = Domain.PiecewiseLinearComplexDomain(vertices=vertices,
                                                  vertexFlags=vertexFlags,
                                                  facets=facets,
                                                  facetFlags=facetFlags,
                                                  holes=holes,
                                                  regions=regions,
                                                  regionFlags=regionFlags)

domain.writePoly("hollow_cylinder")  # Generates hollow_cylinder.poly

analyticalSolution = None

viscosity     = 8.9e-4  #kg/(m*s)
density       = 998.2   #kg/m^3
gravity       = 9.8     #m/s^2
beta          = 0.0 #density*gravity*4.524e-10
m_per_s_by_m_per_d = 1.1574074e-5
lengthScale   = 1.0     #m
timeScale     = 1.0     #d #1.0/sqrt(g*lengthScale)
#make non-dimensional

T = 5e-4       # m²/s
b = height     # 3.0 m
Ks = T / b     # m/s

# Storage Zone
permeability1  = Ks *viscosity/(gravity*density)  #m^2
#permeability1  = (0.00504)*viscosity/(gravity*density)  #m^2
thetaS1        = 0.4   #-
thetaR1        = 0.05   #-
mvg_alpha1     = 8   #1/m
mvg_n1         = 2.4
mvg_m1         = 1.0 - 1.0/mvg_n1
dimensionless_conductivity1  = (timeScale*density*gravity*permeability1/(viscosity*lengthScale))

#print 'Ks',dimensionless_conductivity
dimensionless_density  = 1.0
dimensionless_gravity  = numpy.array([0.0,
                                        -1.0,
                                        0.0])
#dimensionless_alpha    = mvg_alpha*lengthScale
nMediaTypes  = 1
alphaVGtypes = numpy.zeros((nMediaTypes+1,),'d')
nVGtypes     = numpy.zeros((nMediaTypes+1,),'d')
thetaStypes  = numpy.zeros((nMediaTypes+1,),'d')
thetaRtypes  = numpy.zeros((nMediaTypes+1,),'d')
thetaSRtypes = numpy.zeros((nMediaTypes+1,),'d')
KsTypes      = numpy.zeros((nMediaTypes+1,2),'d')

for i in range(nMediaTypes+1):
    if i==1:
        alphaVGtypes[i] = mvg_alpha1
        nVGtypes[i]     = mvg_n1
        thetaStypes[i]  = thetaS1
        thetaRtypes[i]  = thetaR1
        thetaSRtypes[i] = thetaStypes[i] - thetaRtypes[i]
        KsTypes[i,:]    = [dimensionless_conductivity1,dimensionless_conductivity1]#m/d?
    else:
        alphaVGtypes[i] = mvg_alpha1
        nVGtypes[i]     = mvg_n1
        thetaStypes[i]  = thetaS1
        thetaRtypes[i]  = thetaR1
        thetaSRtypes[i] = thetaStypes[i] - thetaRtypes[i]
        KsTypes[i,:]    = [dimensionless_conductivity1,dimensionless_conductivity1]#m/d?


galerkin = False
#useSeepageFace = True

#def getSeepageFace(flag):
#    if useSeepageFace:
#        if flag == boundaryTags['drain']:
#            return 1
#        else:
#            return 0
#    else:
#        return 0




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
                                     diagonal_conductivity=True,
                                     STABILIZATION_TYPE=2,
                                     ENTROPY_TYPE=1,
                                     LUMPED_MASS_MATRIX= False ,
                                     FCT=False, #False ,#True,
                                     num_fct_iter=0,
                                     # FOR ENTROPY VISCOSITY
                                     cE=1.0,
                                     uL=0.0,
                                     uR=1.0,
                                     # FOR ARTIFICIAL COMPRESSION
                                     cK=1.0,
                                     # OUTPUT quantDOFs
                                     outputQuantDOFs=False,
                                     storavity= 1e-4) #,
#                                     getSeepageFace=getSeepageFace)
#galerkin = False

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

#G= [1.8, 1.2,1]

galerkin = False

######THEIS Head at pump################
from scipy.special import exp1
import numpy as np

# Parameters
Q = 1e-4        # pumping rate [m^3/s]
T = 5e-4        # transmissivity [m^2/s]
S = 1e-3        # storativity [-]
h0 = 3.0        # initial head [m]
#pump_radius = 0.6  # distance from center to pump wall

def theis_head_at_pump(t):
    if t <= 0.0:
        return h0
    u = (pump_radius**2 * S) / (4 * T * t)
    drawdown = (Q / (4 * np.pi * T)) * exp1(u)
    return h0 - drawdown

#G=[300.0,40.0,1.0]

#pondingPressure= 0.5
Water_Table =40.0

def getDBC_3D_TheisPumpHead(x, flag):
    if flag == boundaryTags['pump']:
        return lambda x, t: theis_head_at_pump(t) - x[2]
    

dirichletConditions = {0:getDBC_3D_TheisPumpHead}
   
    


# def getDBC_2D_Richards_Shock(x,flag):
# #    return None
#     if 13.99 < x[0]< 16.001:
#         if x[1]>= 20.0 and x[1]<= 39.0:
#             return lambda x,t: -2e-5
#     if (x[0] == 0.0 or
#         x[0] == L[0]):
#         return lambda x,t: (x[1] - L[1]) *dimensionless_gravity[1]*dimensionless_density
#dirichletConditions = {0:getDBC_2D_Richards_Shock}

h0 = 3.0
class ShockIC_2D_Richards:
    def uOfXT(self, x, t):
        bc=getDBC_3D_TheisPumpHead(x,0)
        if bc != None:
            return bc(x,t)
        else:
            return h0 - x[2]
    
    # def uOfXT(self,x,t):
    #     bc=getDBC_2D_Richards_Shock(x,0)
    #     if bc != None:
    #         return bc(x,t)
    #     else:
    #         return (x[1] - L[1]) * dimensionless_gravity[1] * dimensionless_density #G[1] - x[1]
           # z = x[1]
           # if z <= Water_T able:
           #     return Water_Table - z  # Pressure head in meters
           # else:
           #     return 0.0  # Some negative value indicating unsaturated zone

initialConditions  = {0:ShockIC_2D_Richards()}


Q = 1e-4          # m³/s (total pumping rate)

A = 2 * np.pi * pump_radius * height
q_flux = -Q / A   # negative for extraction

q_inj = 2e-5  # example injection rate in m/s or appropriate units
def getFBC_2D_Richards_Shock(x,flag):
    if flag == boundaryTags['pump']:
        return lambda x, t: q_flux
    else:
        return None
    
#fluxBoundaryConditions = {0:'noFlow'}
fluxBoundaryConditions = {0:getFBC_2D_Richards_Shock}



#        return lambda x,t:0.0

advectiveFluxBoundaryConditions =  {0:getFBC_2D_Richards_Shock}

diffusiveFluxBoundaryConditions = {0:{0:getFBC_2D_Richards_Shock}}

T = 0.2/timeScale

