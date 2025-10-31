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
from domain_rg import *


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

# # Finite element sapce #
# pDegree_tadr=1
# useBernstein=False
# useHex=False

# # quadrature order #
# tadr_quad_order = 2*pDegree_tadr+1
# #from gw_cl_domain import *
# SSPOrder = 2

class MyCoefficients(TADR.Coefficients):
    def attachModels(self, modelList):
        """
        Attach the Richards model to TADR and pass the velocity field (grad(u)).
        """
        self.model = modelList[self.modelIndex]
        self.vModel = modelList[self.V_model]
        # Attach Richards model
        logEvent(f"[TADR.attachModels] Bound TADR model at idx={self.modelIndex}, "
                 f"type={type(self.model).__name__}, id={id(self.model)}")
        logEvent(f"[TADR.attachModels] Coeff==model.coeff? {self.model.coefficients is self}")

#        flowModel = modelList[0]  # Richards model

        for i, model in enumerate(modelList):
            logEvent(f"Model {i}: {model.name}, Type: {type(model)}")
        self.q_v    = self.vModel.q['velocity']
        self.ebqe_v = self.vModel.ebqe['velocity']

        # self.q_v = flowModel.q['velocity']
        # self.ebqe_v = flowModel.ebqe['velocity']
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
alpha_L= 0.3
alpha_T= 0.1*alpha_L
Dm= 0.01 #18.86e-6
coefficients = MyCoefficients(
    alpha_L=alpha_L,
    alpha_T=alpha_T,
    Dm=Dm,
    rho_fw=1000.0,
    rho_sw=1025.0,
    V_model=0,  
    checkMass=checkMass,
    FCT=True,
    LUMPED_MASS_MATRIX=True, 
    STABILIZATION_TYPE=2,
    diagonal_conductivity= True, 
    ENTROPY_TYPE='LOG', 
    cE=0.1, cK=1.0, physicalDiffusion=0.0) #ct.physicalDiffusion) 
coefficients.variableNames=['u']

#####################
# INITIAL CONDITION #
#####################
###########################
# Defining Initial Conditions Functions
# Initially will start at all freshwater
# and then go to a steady state
###########################

# ---- geometry extents reused from flow ----


cin = 3.0  # inlet concentration during pulses

# ---- pulse schedule (in DAYS) ----
# 30 minutes = 0.5/24 days
pulse_len = 2.0/24.0
# three pulses starting at day 0, day 1, day 2 (adjust as you like)
pulse_starts = [0.0, 10.0, 40.0]

def in_pulse(t):
    """Return True if time t (in days) is inside any pulse window."""
    for t0 in pulse_starts:
        if t0 <= t <= t0 + pulse_len:
            return True
    return False

# -------------------------
# Initial condition (clean)
# -------------------------
class ConstantIC:
    def __init__(self, cval=0.0):
        self.cval = cval
    def uOfXT(self, x, t):
        # If a pulse is active and we're on the top-middle strip, start at cin
        if (in_pulse(t) and
            x[1] >= 1.7 and (G[0]/3.0 <= x[0] <= 2.0*G[0]/3.0)):
            return cin
        return self.cval
    def uOfX(self, x):
        return self.cval

initialConditions = {0: ConstantIC(0.0)}

# -------------------------------------------
# Dirichlet BC: only on ponded top sub-segment,
# and only during pulse windows
# -------------------------------------------
def getConcDirichletBC(x, flag):
    if flag == boundaryTags['top'] and (G[0]/3.0 <= x[0] <= 2.0*G[0]/3.0):
        return lambda x, t: (cin if in_pulse(t) else 0.0)
    return None

dirichletConditions = {0: getConcDirichletBC}

# ---------------------------------------------------------
# Flux BCs:
#  - Advective: use velocity ('setFlow')
#  - Diffusive: zero everywhere except where Dirichlet applies
# ---------------------------------------------------------
fluxBoundaryConditions = {0: 'setFlow'}   # let velocity drive inflow/outflow

# leave advectiveFluxBoundaryConditions unset (or empty), do NOT set {0: None}
# advectiveFluxBoundaryConditions = {}

def getZeroDiffusiveFlux(x, flag):
    # Where Dirichlet may apply (the inlet strip), return None so we don't overconstrain
    if flag == boundaryTags['top'] and (G[0]/3.0 <= x[0] <= 2.0*G[0]/3.0):
        return None
    return lambda x, t: 0.0
advectiveFluxBoundaryConditions = {0: getZeroDiffusiveFlux}  # use 'setFlow' behavior
diffusiveFluxBoundaryConditions = {0: {0: getZeroDiffusiveFlux}}