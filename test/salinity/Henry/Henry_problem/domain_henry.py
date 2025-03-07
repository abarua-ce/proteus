from proteus import *
import numpy
from math import *
from proteus import Domain
from proteus import Norms
from proteus import Profiling 
from proteus import Context 
#from proteus.mprans import TADR
import numpy as np
import math 
from proteus.Profiling import logEvent
from proteus.default_p import *
from proteus.richards import Richards

details = """
The Dispersive Henry Benchmark problem (Arbarca et.al. 2007) and the modified Henry Benchmark
Problem of Simpson and Clement (2004):
Variable Density Groundwater Flow and Transport
in homogeneous porous medium
The primitive variables are the mass fraction (mf) 
and the freshwater pressure head (h)
"""

###########################
# Defining Problem Flags, Spatial, Temporal Domain Data, Transport Model
#  Numerical Information common to both nfiles
###########################

# Flags that set the variable viscosity, consistent velocity algorithm
# and velocity postprocessing
#variable_viscosity = 0 # 0 - None, 1 - Variable Viscosity
#consistentVelocityFlag=0 # 0, No Consistent Velocity Formulation, 1 Consistent
#unstructure=True
from proteus import MeshTools
partitioningType = MeshTools.MeshParallelPartitioningTypes.node

#spatial domain Omega = [0,2.0] x [0,1.0]
L=(2.0,1.0) # m, Henry
nd = 2
nnx= 80
nny=40

domain = Domain.RectangularDomain(L=[L[0],L[1]],name="gwdtrans_domain",units="m")
polyfile="gwvdtransport_domain_2d"
domain.writePoly(polyfile)
triangleOptions = "q30Dena0.00001"
domain.MeshOptions.nnx = nnx
domain.MeshOptions.nny = nny
domain.MeshOptions.triangleFlag=0

#else:
#    nnx = 80
#    nny = 40
 

#T =36000.0 #100.0 #12000 # 6.0e2*6..,1.0e4 #time scale, s
#nDTout =10 #100  # output timesteps to use
#DT=1

