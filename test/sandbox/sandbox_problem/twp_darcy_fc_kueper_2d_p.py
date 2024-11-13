from proteus import *
from proteus.default_p import *
#from proteus.TransportCoefficients import *
from proteus.TwophaseDarcyCoefficients import *

#from proteus.SubsurfaceTransportCoefficients import *

from impes_kueper_setParams import *

analyticalSolutions = None

g=[0, -9.8 , 0]
# Define sparseDiffusionTensors dynamically
def buildSparseDiffusionTensors(is_diagonal, Ksw_block, nd):
    sparseDiffusionTensors = {}
    for i, Ksw in enumerate(Ksw_block):
        if is_diagonal:
            sparseDiffusionTensors[(i, 0, 0)] = (
                np.arange(nd + 1, dtype='i'), np.arange(nd, dtype='i')
            )
        else:
            sparseDiffusionTensors[(i, 0, 0)] = (
                np.arange(nd**2 + 1, step=nd, dtype='i'), 
                np.array([list(range(nd)) for _ in range(nd)], dtype='i')
            )
    return sparseDiffusionTensors

# Create the sparseDiffusionTensors dictionary
sparseDiffusionTensors = buildSparseDiffusionTensors(is_diagonal=True, Ksw_block=Ksw_block, nd=nd)


coefficients = TwophaseDarcy_fc(
    Ksw=Ksw,                      # Imported Ksw
    rhon=rhon,                    # Imported density of non-wetting phase
    rhow=rhow,                    # Imported density of wetting phase
    g=g,
    vg_alpha=mvg_alpha,            # Imported van Genuchten alpha parameter
    bc_lambda=bc_lambda,          # Imported Brooks-Corey lambda
    bc_pd=bc_pd,                  # Imported Brooks-Corey entry pressure
    omega=omega,                  # Imported porosity
    Sw_max=sw_max_block[0],                # Imported maximum saturation
    Sw_min=sw_min_block[0],                # Imported residual saturation
    #density_w_parameters=density_w_parameters,  # Imported density parameters for wetting phase
    #density_n_parameters=density_n_parameters,  # Imported density parameters for non-wetting phase
    sparseDiffusionTensors=sparseDiffusionTensors  # Set based on whether diffusion is isotropic or anisotropic
)


# coefficients = TwophaseDarcy_fc(Ksw=Ksw,
#                                    rhon=rhon,
#                                    rhow=rhow,
#                                    g=g,
#                                    vg_alpha=5.0,
#                                    bc_lambda=bc_lambda,
#                                    bc_pd = bc_pd)#,
#                                    mvg_n = mvg_n,
#                                    mvg_m = mvg_m,
#                                    omega=omega,
#                                    mun=mun,
#                                    muw=muw)#,
#                                   model=model,
#                                   setParamsFunc=setParams)

# if useHet:
#     coefficients = TwophaseFFDarcyFCHet(Ksw=Ksw,
#                                         rhon=rhon,
#                                         rhow=rhow,
#                                         g=g,
#                                         mvg_alpha=mvg_alpha,
#                                         bc_lambda=bc_lambda,
#                                         bc_pd = bc_pd,
#                                         mvg_n = mvg_n,
#                                         mvg_m = mvg_m,
#                                         omega=omega,
#                                         mun=mun,
#                                         muw=muw,
#                                         model=model,
#                                         setParamsFunc=setParams)
# else:
#     coefficients = TwophaseFFDarcyFC(Ksw=Ksw,
#                                    rhon=rhon,
#                                    rhow=rhow,
#                                    g=g,
#                                    mvg_alpha=mvg_alpha,
#                                    bc_lambda=bc_lambda,
#                                    bc_pd = bc_pd,
#                                    mvg_n = mvg_n,
#                                    mvg_m = mvg_m,
#                                    omega=omega,
#                                    mun=mun,
#                                    muw=muw,
#                                    model=model)
    
# Define any undefined parameters or flags here
top = 1.0          # Set the y-coordinate for the top boundary
right = 1.0        # Set the x-coordinate for the right boundary
slit_is_top = True # Boolean flag indicating if the slit is at the top
FudgeFactor = 0.05 # Adjust this value as required
open_bottom = True # Boolean flag indicating if the bottom boundary is open

# Now include the Dirichlet boundary conditions functions
def getDBC_sw(x):
    # constant saturation over slit
    if x[1] == top:
        if slit_is_top or (x[0] >= right / 3.0 and x[0] <= 2.0 * right / 3.0):
            return lambda x, t: 1.0 - FudgeFactor
    if open_bottom and x[1] == 0.0:
        return lambda x, t: 1.0 - FudgeFactor
    elif x[0] in [0.0, right]:  # open sides
        return lambda x, t: FudgeFactor

def getDBC_psiw(x):
    # constant head over slit
    if x[1] == top:
        if slit_is_top or (x[0] >= right / 3.0 and x[0] <= 2.0 * right / 3.0):
            return lambda x, t: 0.1
    if open_bottom and x[1] == 0.0:
        return lambda x, t: 0.0
    elif x[0] in [0.0, right]:  # open sides
        return lambda x, t: -x[1]  # 0.0

# Dirichlet and initial conditions setup
dirichletConditions = {0: getDBC_sw, 1: getDBC_psiw}

class sw_IC:
    def uOfXT(self, x, t):
        if x[1] == top and (slit_is_top or (x[0] >= right / 3.0 and x[0] <= 2.0 * right / 3.0)):
            return 1.0 - FudgeFactor
        else:
            return FudgeFactor

class psiw_IC:
    def uOfXT(self, x, t):
        if x[1] == top and (slit_is_top or (x[0] >= right / 3.0 and x[0] <= 2.0 * right / 3.0)):
            return 0.1
        else:
            return 0.0

initialConditions = {0: sw_IC(), 1: psiw_IC()}

# Define boundary conditions for fluxes
fluxBoundaryConditions = {0: 'noFlow', 1: 'noFlow'}

# Commented sections for custom advective and diffusive flux boundary conditions
advectiveFluxBoundaryConditions = {}
diffusiveFluxBoundaryConditions = {0: {}, 1: {}}


# #now define the Dirichlet boundary conditions
# def getDBC_sw(x):
#     #constant saturation over slit
#     if x[1] == top:
#         if (slit_is_top or
#             (x[0] >= right/3.0 and
#              x[0] <= 2.0*right/3.0)):
#             return lambda x,t: 1.0-FudgeFactor
#     if open_bottom:
#         if x[1] == 0.0:
#             return lambda x,t: 1.0-FudgeFactor
#     else:#open sides
#         if x[0] in [0.0,right]:
#             return lambda x,t: FudgeFactor

# def getDBC_psiw(x):
#     #constant head over slit
#     if x[1] == top:
#         if (slit_is_top or
#             (x[0] >=right/3.0 and
#             x[0] <= 2.0*right/3.0)):
#             return lambda x,t: 0.1
#     if open_bottom:
#         if x[1] == 0.0:
#             return lambda x,t: 0.0
#     else:#open sides
#         if x[0] in [0.0,right]:
#             return lambda x,t: -x[1]#0.0

# dirichletConditions = {0:getDBC_sw,1:getDBC_psiw}

# class sw_IC:
#     def __init__(self):
#         pass
#     def uOfXT(self,x,t): 
#         if (x[1] == top and
#             (slit_is_top or
#             (x[0] >= right/3.0 and
#              x[0] <= 2.0*right/3.0))):
#             return 1.0-FudgeFactor
#         else:
#             return FudgeFactor

# class psiw_IC:
#     def __init__(self):
#         pass
#     def uOfXT(self,x,t):
#         if (x[1] == top and
#             (slit_is_top or
#             (x[0] >=right/3.0 and
#              x[0] <= 2.0*right/3.0))):
#             return 0.1
#         else:
#             return 0.0
	
# initialConditions  = {0:sw_IC(),1:psiw_IC()}

# fluxBoundaryConditions = {0:'noFlow',1:'noFlow'}


# def get_w_AFBC(x):
#     #no flow outside  of slit and on bottom
#     if x[1] == top:
#         if (not slit_is_top or
#             (x[0] < right/3.0 or
#             x[0] > 2.0*right/3.0)):
#             return lambda x,t: 0.0
#     if x[1] == 0.0:
#         return lambda x,t: 0.0

# def get_n_AFBC(x):
#     #no flow outside  of slit and on bottom
#     if x[1] == top:
#         if (not slit_is_top or
#             (x[0] < right/3.0 or
#             x[0] > 2.0*right/3.0)):
#             return lambda x,t: 0.0
#     if x[1] == 0.0:
#         return lambda x,t: 0.0

# advectiveFluxBoundaryConditions =  {0:get_w_AFBC,
#                                     1:get_n_AFBC}
#advectiveFluxBoundaryConditions =  {}
#diffusiveFluxBoundaryConditions = {0:{},1:{}}


