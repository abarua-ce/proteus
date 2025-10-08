from __future__ import absolute_import
from builtins import object
import numpy as np
from proteus import *
from proteus.default_p import *
from proteus.ctransportCoefficients import smoothedHeaviside

from proteus import Domain
from proteus import Norms
from proteus import Profiling 
from proteus import Context 
from proteus.mprans import TADR

from proteus.Profiling import logEvent

from math import *
try:
    from .domain_henry import *
except:
    from domain_henry import *

#T =4800.0 #100.0 #12000 # 6.0e2*6..,1.0e4 #time scale, s
#nDTout =10 #0  # output timesteps to use

################################## Inputs ########################################

physicalDiffusion = 0.0 #1.0e-8 # m2/s, molecular diffusion coefficient
refinement = 0
#unstructured = False



# SHOCK CAPTURING PARAMETERS #
shockCapturingFactor_tadr=0.2
lag_shockCapturing_tadr=True

# number of space dimensions #
#nd=ct.nd

# General parameters #
parallel = False # if True use PETSc solvers
linearSmoother = None
checkMass = False

# Finite element sapce #
pDegree_tadr=1
useBernstein=False
useHex=False

# quadrature order #
tadr_quad_order = 2*pDegree_tadr+1
#from gw_cl_domain import *
SSPOrder = 2

class MyCoefficients(TADR.Coefficients):
    def attachModels(self, modelList):
        """
        Attach the Richards model to TADR and pass the velocity field (grad(u)).
        """
        self.model = modelList[self.modelIndex]
        # Attach Richards model
        logEvent(f"[TADR.attachModels] Bound TADR model at idx={self.modelIndex}, "
                 f"type={type(self.model).__name__}, id={id(self.model)}")
        logEvent(f"[TADR.attachModels] Coeff==model.coeff? {self.model.coefficients is self}")

        flowModel = modelList[0]  # Richards model
        
        for i, model in enumerate(modelList):
            logEvent(f"Model {i}: {model.name}, Type: {type(model)}")
        self.q_v = flowModel.q['velocity']
        self.ebqe_v = flowModel.ebqe['velocity']
        # Log for debugging
        logEvent(f"Richards velocity (q_v): mean={self.q_v.mean()}, min={self.q_v.min()}, max={self.q_v.max()}")
        logEvent(f"Richards boundary velocity (ebqe_v): mean={self.ebqe_v.mean()}, min={self.ebqe_v.min()}, max={self.ebqe_v.max()}")
        logEvent("Attached Richards model to TADR coefficients.")


LevelModelType = TADR.LevelModel
logEvent = Profiling.logEvent
#name=soname
a0= 18.8571e-6 #e-6 
def a(x):
    return np.array([[a0,0.0],[0.0,a0]])
aOfX = {0:a}
alpha_L= 0.1
alpha_T= 0.1*alpha_L
Dm= 18.86e-6
coefficients = MyCoefficients(
    alpha_L=alpha_L,
    alpha_T=alpha_T,
    Dm=Dm,
    V_model=0,  
    checkMass=checkMass,
    FCT=True,
    LUMPED_MASS_MATRIX=True, 
    STABILIZATION_TYPE=2,
    diagonal_conductivity= True, 
    ENTROPY_TYPE='LOG', 
    cE=0.1, cK=1.0, physicalDiffusion=0.0) #ct.physicalDiffusion) 
coefficients.variableNames=['u']

# Density Relation
rho_fw =  1000.0 # density of freshwater
rho_sw = 1025.0 # 1026.0 density of saltwater
mf_fw = 0.0 # saline mass fraction, freshwater
mf_sw = 10.0# .0357 # saline mass fraction, saltwater
epsilon = (rho_sw - rho_fw)/rho_fw 
drho_dmf =epsilon*rho_fw #700. #epsilon*rho_fw # density change with mass fraction

# Dynamic Viscosity Relation
nu = 1.0e-3

#####################
# INITIAL CONDITION #
#####################
###########################
# Defining Initial Conditions Functions
# Initially will start at all freshwater
# and then go to a steady state
###########################


class constantIC:
    def __init__(self,cval=0.0):
        self.cval = cval
    def uOfXT(self,x,t):
        return self.cval
    def uOfX(self,x):
        return self.cval
#        return self.cval


def getHenryConcDirichletBCs(x,tag):
    if x[0]<1.e-8:     
        return lambda x,t: mf_fw 
    elif abs(x[0]-L[0])<1.e-8:
        return lambda x,t: mf_sw

# Henry Problem Mass Flux Boundary Conditions

def getzeroMassDiffusiveFluxBCs(x,tag):
    if ((x[0] > 1.e-8 and abs(x[0] -L[0])> 1.e-8) and 
        (x[1] < 1.0e-8 or  abs(x[1]- L[1]) < 1.0e-8) ):
        return lambda x,t: 0.0
    else: 
        pass

def getzeroBCs(x,tag):
    return lambda x,t: 0.0
   
def getzeroMassAdvectiveFluxBCs(x,tag):
    if ((x[0] > 1.0e-8 and x[0] < L[0]-1.0e-8) and 
       ( x[1] < 1.0e-8 or x[1] > L[1]-1.0e-8)):
        return lambda x,t: 0.0
    

initialConditions ={0:constantIC(mf_fw)} 
dirichletConditions ={0:getHenryConcDirichletBCs} 
fluxBoundaryConditions = {0:'setFlow'}
advectiveFluxBoundaryConditions ={0:getzeroMassDiffusiveFluxBCs}
diffusiveFluxBoundaryConditions ={0:{0:getzeroMassDiffusiveFluxBCs}}#{0:{0:getzeroBCs}}#{0:{0:getzeroMassDiffusiveFluxBCs}}

