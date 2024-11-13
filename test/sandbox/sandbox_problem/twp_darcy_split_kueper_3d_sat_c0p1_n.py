from proteus import *
from proteus.default_n import *
from twp_darcy_split_kueper_3d_sat_p import *

#general type of integration (Forward or to SteadyState)
timeIntegrator = ForwardIntegrator
#type of time integration formula
#timeIntegration = BackwardEuler
# stepController = FixedStep
#DT=1.0e1 
#stepController = Min_dt_controller
timeIntegration = FLCBDF
stepController = FLCBDF_controller
rtol_u[1] = 1.0e-3
rtol_u[2] = 1.0e-3
atol_u[1] = 1.0e-3
atol_u[2] = 1.0e-3
#runCFL = 1000.0
#runCFL = 200.0
runCFL=None
DT = None
nDTout = 50#int(T/DT)
print "nDTout",nDTout
femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,3)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

#elementQuadrature = SimplexLobattoQuadrature(nd,1)
#elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,1)

nn=3
nLevels = 1
triangleOptions = "pAfena%e" % ((top*right*top*0.1/6.0)*0.001,)



subgridError = None
subgridError = AdvectionDiffusionReaction_ASGS(coefficients,nd,stabFlag='2',lag=True)


massLumping=False

#shockCapturing = None
shockCapturing = ResGrad_SC(coefficients,nd,shockCapturingFactor=0.5,lag=True)

multilevelNonlinearSolver  = Newton

levelNonlinearSolver = Newton
maxNonlinearIts = 20
maxLineSearches = 10

fullNewtonFlag = True

tolFac = 0.001

nl_atol_res = 1.0e-7

matrix = SparseMatrix

multilevelLinearSolver = PETSc#PETSc LU

levelLinearSolver = PETSc#PETSc LU
#pick number of layers to use in overlap 
#"-ksp_type cg -pc_type asm -pc_asm_type basic -ksp_atol  1.0e-10 -ksp_rtol 1.0e-10 -ksp_monitor_draw" or
#-pc_type lu -pc_factor_mat_solver_package
nLayersOfOverlapForParallel = 1
#type of partition
parallelPartitioningType = MeshParallelPartitioningTypes.node
#parallelPartitioningType = MeshParallelPartitioningTypes.element
numericalFluxType = Advection_DiagonalUpwind_Diffusion_IIPG_exterior


linTolFac = 0.0001

conservativeFlux = None
