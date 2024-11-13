from __future__ import print_function
from __future__ import absolute_import
from builtins import range
from proteus import *
from proteus.default_n import *
from tadr_CCS_p import *

multilevelNonlinearSolver  = Newton
if ct.STABILIZATION_TYPE in ['Galerkin', 'VMS']: #SUPG
    levelNonlinearSolver = Newton

# SHOCK CAPTURING PARAMETERS #
shockCapturingFactor_tadr=0.2
lag_shockCapturing_tadr=True



# # General parameters #
# parallel = False # if True use PETSc solvers
# linearSmoother = None
# checkMass = False

# Finite element sapce #
pDegree_tadr=1
useBernstein=False
useHex=False

# quadrature order #
tadr_quad_order = 2*pDegree_tadr+1

# parallel partitioning info #
from proteus import MeshTools
partitioningType = MeshTools.MeshParallelPartitioningTypes.node
# create mesh #
nn=nnx=(2**ct.refinement)*10+1
nny=1
if nd == 2:
    nny=(nnx-1)//10+1 


domain.MeshOptions.nn = nn
domain.MeshOptions.nnx = nnx
domain.MeshOptions.nny = nny
domain.MeshOptions.nnz = nnz
domain.MeshOptions.triangleFlag=0


soname="tadr_level_"+repr(ct.refinement)


multilevelNonlinearSolver  = Newton
if ct.STABILIZATION_TYPE in ['Galerkin', 'VMS']: #SUPG
    levelNonlinearSolver = Newton
    maxLineSearches = 0
    fullNewtonFlag = True
    updateJacobian = True
    timeIntegration = BackwardEuler_cfl
elif ct.STABILIZATION_TYPE=='TaylorGalerkinEV':
    levelNonlinearSolver = TwoStageNewton
    fullNewtonFlag = False
    updateJacobian = False
    timeIntegration = BackwardEuler_cfl
else:
    fullNewtonFlag = False
    updateJacobian = False
    timeIntegration = TADR.RKEV # SSP33
    if ct.LUMPED_MASS_MATRIX==True:
        levelNonlinearSolver = ExplicitLumpedMassMatrix
    else:
        levelNonlinearSolver = ExplicitConsistentMassMatrixForVOF

stepController = Min_dt_controller
runCFL = ct.cfl
timeOrder = ct.SSPOrder
nStagesTime = ct.SSPOrder

if useHex:
    hex=True
    quad=True
    if pDegree_tadr == 1:
        femSpaces = {0:C0_AffineLinearOnCubeWithNodalBasis}
    elif pDegree_tadr == 2:
        if useBernstein:
            femSpaces = {0:C0_AffineBernsteinOnCube}
        else:
            femSpaces = {0:C0_AffineLagrangeOnCubeWithNodalBasis}
    else:
        print("pDegree = %s not recognized " % pDegree_tadr)
    elementQuadrature = CubeGaussQuadrature(nd,tadr_quad_order)
    elementBoundaryQuadrature = CubeGaussQuadrature(nd-1,tadr_quad_order)
else:
    if pDegree_tadr == 1:
        femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}
    elif pDegree_tadr == 2:
        femSpaces = {0:C0_AffineQuadraticOnSimplexWithNodalBasis}
    else:
        print("pDegree = %s not recognized " % pDegree_tadr)
    elementQuadrature = SimplexGaussQuadrature(nd,tadr_quad_order)
    elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,tadr_quad_order)

numericalFluxType = TADR.NumericalFlux
#numericalFluxType = DoNothing
#numericalFluxType = Advection_DiagonalUpwind_IIPG_exterior # PERIODIC

shockCapturing = TADR.ShockCapturing(coefficients,nd,shockCapturingFactor=shockCapturingFactor_tadr,lag=lag_shockCapturing_tadr)

matrix = SparseMatrix
if parallel:
    multilevelLinearSolver = KSP_petsc4py
    levelLinearSolver = KSP_petsc4py
    linear_solver_options_prefix = 'tadr_'
    linearSolverConvergenceTest = 'r-true'
else:
    multilevelLinearSolver = LU
    levelLinearSolver = LU

if checkMass:
    auxiliaryVariables = [MassOverRegion()]

tnList=[0.,1E-6]+[float(n)*ct.T/float(ct.nDTout) for n in range(1,ct.nDTout+1)]
