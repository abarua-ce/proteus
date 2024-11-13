#from proteus import Domain
#from proteus import Norms
#from proteus import Profiling 
#from proteus import Context 

from __future__ import absolute_import
from builtins import object
import numpy as np
from proteus import *
from proteus.mprans import TADR
from proteus.default_p import *
from proteus.ctransportCoefficients import smoothedHeaviside
import numpy as np
import math 

ct=Context.Options([
    # General parameters #
    ("physicalDiffusion", 0.0, "isotropic diffusion coefficient"),
    #("problem",0,"0: 1D problem with periodic BCs. 1: 2D rotation of zalesak disk"),
    ("nd",2,"space dimension"),
    ("T",0.1,"Final time"),
    ("nDTout",1,"Number of time steps to archive"),
    ("refinement",0,"Level of refinement"),
    ("unstructured",False,"Use unstructured mesh. Set to false for periodic BCs"),
    # Choice of numerical method #
    ("STABILIZATION_TYPE","Galerkin","Galerkin, VMS, TaylorGalerkinEV, EntropyViscosity, SmoothnessIndicator, Kuzmin, "),
    ("LUMPED_MASS_MATRIX",False,"Flag to lumped the mass matrix"),
    ("ENTROPY_TYPE","POWER","POWER, LOG"),
    ("FCT",True,"Use Flux Corrected Transport"),
    # Numerical parameters #
    ("cE",0.1,"Entropy viscosity constant"),
    ("cK",1.0,"Artificial compression constant"),
    ("cfl",0.5,"Target cfl"),
    ("SSPOrder",2,"SSP method of order 1, 2 or 3")
],mutable=True)


LevelModelType = TADR.LevelModel
logEvent = Profiling.logEvent
soname="tadr_level_"+repr(ct.refinement)
name= soname
a0= 0.02
def a(x):
    return np.array([[a0,0.0],[0.0,a0]])
aOfX = {0:a}

nd= 2

domain = Domain.PlanarStraightLineGraphDomain()

boundaries = [f'marker_{i}' for i in range(1, 17)]

boundaryTags=dict([(key,i+1) for (i,key) in enumerate(boundaries)])

# Read the values from the text file
file_path = 'points_array_try.txt'  # Replace with your actual file path
with open(file_path, 'r') as file:
    lines = file.readlines()

# Extract the values for the first and third columns
selected_values = [[float(line.split(',')[0]), float(line.split(',')[2])] for line in lines]

# Create a NumPy array from the selected values
filtered_points = np.array(selected_values)

#vertices= filtered_points.tolist()

# Rescale to [-1, 1] for both x and y
min_vals = filtered_points.min(axis=0)
max_vals = filtered_points.max(axis=0)
scaled_points = 2 * (filtered_points - min_vals) / (max_vals - min_vals) - 1
vertices = scaled_points.tolist()



vertexFlags = []
marker= np.loadtxt('vertexflags.txt', dtype=int).tolist()
for value in marker:
    formatted_value = f"marker_{value}"
    vertexFlags.append(boundaryTags[formatted_value])

segments=np.loadtxt('segments.txt', dtype=int).tolist()

segmentFlags = []
flags_segments= np.loadtxt('segmentsFlags.txt', dtype=int).tolist()

for value in flags_segments:
    marker_number= int(value)
    formatted_value = f"marker_{marker_number}"
    segmentFlags.append(boundaryTags[formatted_value])

#holes= [[1.5,1.5]]
#regions=[[0.1,0.1], [1.5, 2.5], [1.5, 4.5], [1.5,0.5]]
#regionFlags=[1,2,3,4]
#regions=[[0.1,0.1], [1.5, 2.5], [1.5, 4.5]]#, [1.5,0.5]]
regions= [[0.8,1.4], #1
          [0.8,1.15], #2
           [0.6,1.0],[1.2,0.9],[1.7,1.0], #3
           [0.4,0.9],[1.1,0.8],[2.2,0.95], #4
           [0.5,0.85],[1.2,0.75],[1.7,0.9], #5
           [0.42,0.75],[1.20,0.65],[1.62,0.75], #6
           [2.61,0.70], #7
           [0.18,0.55],[0.60,0.60], #8
           [0.30,0.65], #9
           [0.21,0.60],[0.60,0.65], #10
           [0.15,0.35],[0.60,0.45], #11
           [0.15,0.20],[0.81,0.25], #12
           [1.83,0.05], #13
           [0.36,0.25], #14
           [0.96,0.90], #15
           [1.44,0.85], #16 
           ]
scaled_regions = 2 * (regions - min_vals) / (max_vals - min_vals) - 1
regions= scaled_regions.tolist()
#regionFlags=[1,2,3]
regionFlags=[1,
            2,
             3,3,3,
             4,4,4,
             5,5,5,
             6,6,6,
             7,
             8,8,
             9,
             10,10,
             11,11,
             12,12,
             13,
             14,
             15,16]


domain = Domain.PlanarStraightLineGraphDomain(vertices= vertices,
                                              vertexFlags=vertexFlags,
                                              segments=segments,
                                              segmentFlags=segmentFlags,
                                              regions = regions,
                                              regionFlags = regionFlags,)
                                               #holes= holes)

domain.writePoly('CCS')

# General parameters #
parallel = False # if True use PETSc solvers
linearSmoother = None
checkMass = False

class MyCoefficients(TADR.Coefficients):
    def attachModels(self,modelList):
        self.model = modelList[0]
        self.rdModel = self.model
        self.q_v = np.zeros(self.model.q[('grad(u)',0)].shape,'d')
        self.ebqe_v = np.zeros(self.model.ebqe[('grad(u)',0)].shape,'d')
        # Explicitly apply the velocity field if necessary
        self.q_v[..., 0] = velx(self.model.q[('x')], 0)
        self.q_v[..., 1] = vely(self.model.q[('x')], 0)


coefficients = MyCoefficients(
    aOfX,
    checkMass=checkMass,
    FCT=ct.FCT,
    LUMPED_MASS_MATRIX=ct.LUMPED_MASS_MATRIX, 
    STABILIZATION_TYPE=ct.STABILIZATION_TYPE, 
    ENTROPY_TYPE=ct.ENTROPY_TYPE, 
    cE=ct.cE, cK=ct.cK, physicalDiffusion=ct.physicalDiffusion) 
coefficients.variableNames=['u']



LevelModelType = TADR.LevelModel
logEvent = Profiling.logEvent
name=soname
a0= 0.02
def a(x):
    return numpy.array([[a0,0.0],[0.0,a0]])
aOfX = {0:a}

#nd=ct.nd



coefficients = MyCoefficients(
    aOfX,
    checkMass=checkMass,
    FCT=ct.FCT,
    LUMPED_MASS_MATRIX=ct.LUMPED_MASS_MATRIX, 
    STABILIZATION_TYPE=ct.STABILIZATION_TYPE, 
    ENTROPY_TYPE=ct.ENTROPY_TYPE, 
    cE=ct.cE, cK=ct.cK, physicalDiffusion=ct.physicalDiffusion) 
coefficients.variableNames=['u']




##################
# VELOCITY FIELD #
##################
def velx(X,t):
    return -1.0
    #if ct.problem in [0,2]:
    #    return 1.0
    #else:        
    #    return -2*pi*(X[1]-0.5)

def vely(X,t):
    return -1.0
    #if ct.problem in [0,2]:
    #    return 0.0
    #else:
    #    return 2*pi*(X[0]-0.5)

velocityFieldAsFunction={0:velx,
                         1:vely}

#galerkin = False

#G=[2.799999952316284,1.4950000047683716]
G=[0.999999965940203, 1.0000000063790924]

#pondingPressure= 5.0

#####################
# INITIAL CONDITION #
#####################
class init_cond(object):
    def uOfXT(self,x,t):
          if x[1]>=1.00:
              if (x[0]>=G[0]/3.0 and x[0]<=2*G[0]/3.0):  
                return 3.0
          else:
              return 0.0
            
        #   if x[1]>= 1.44:
        #         if (x[0] >= G[0]/3.0 and
        #         x[0] <= 2.0*G[0]/3.0):
        # 	        return lambda x,t: 3.0
        #     else:
            
        #     if x[1] == 0.0:
        #         return lambda x,t: 0.0
        #     if (x[0] == 0.0 or
        #         x[0] == G[0]):
        #         return lambda x,t: 0.0



    
initialConditions  = {0:init_cond()}
analyticalSolution = {0:init_cond()}

#######################
# BOUNDARY CONDITIONS #
#######################
def getDBC(x,flag):
    if x[1]>=1.00:
              if (x[0]>=G[0]/3.0 and x[0]<=2*G[0]/3.0):  
                #print(f"Dirichlet BC applied at x={x} with value=3.0")
                return lambda x,t: 3.0
    else:
              #print(f"Dirichlet BC applied at x={x} with value=0.0")
              return lambda x,t: 0.0
    #return lambda x,t: 0.0
    
dirichletConditions = {0:getDBC}

def zeroadv(x,flag):
    #if flag == 1:  # Apply a non-zero flux on a specific boundary
    #    print(f"Setting non-zero advective flux on boundary {flag}")
    #    return lambda x, t: 1.0  # Example non-zero flux
    #print(f"zeroadv called with x={x}, flag={flag}")
    return None

advectiveFluxBoundaryConditions =  {0:zeroadv}

def zerodiff(x,flag):
    #print(f"zerodiff called with x={x}, flag={flag}")
    return lambda x,t:0.0 # None

fluxBoundaryConditions = {0:'outFlow'}
diffusiveFluxBoundaryConditions = {0:{0:zerodiff}}
#diffusiveFluxBoundaryConditions = {0:{}}


