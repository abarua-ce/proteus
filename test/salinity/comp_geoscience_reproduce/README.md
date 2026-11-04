# Computational Geosciences runs

Richards (flow) STABILIZATION_TYPE 2 + FCT, TADR (transport) STABILIZATION_TYPE 5 + FCT.
The low-order runs (Zalesak FCT=False, MMS MMS_FCT=0, henry/low_order/) turn FCT off in both.
Run each command from the folder named in its heading.

## 1. Zalesak's disk (`zalesak/`)

```bash
# FCT
mpiexec -n 64 parun thelper_tadr_p.py thelper_tadr_n.py -l3 \
  -C "problem=1 nd=2 T=1.0 nDTout=1 refinement=4 STABILIZATION_TYPE=5 cfl=0.5 SSPOrder=2 FCT=True parallel=True" \
  -P "-tadr_ksp_type preonly -tadr_pc_type lu -tadr_pc_factor_mat_solver_type superlu_dist"
# low order
mpiexec -n 64 parun thelper_tadr_p.py thelper_tadr_n.py -l3 \
  -C "problem=1 nd=2 T=1.0 nDTout=1 refinement=4 STABILIZATION_TYPE=5 cfl=0.5 SSPOrder=2 FCT=False parallel=True" \
  -P "-tadr_ksp_type preonly -tadr_pc_type lu -tadr_pc_factor_mat_solver_type superlu_dist"
```

`refinement=5` for the finer mesh.

## 2. Manufactured solution (`mms/`)

```bash
RICHARDS_EDGE_VELOCITY=1 MMS_N=40 MMS_FCT=1 mpiexec -np 1 parun mms_so.py -l3 \
  -P "-richards_ksp_type preonly -richards_pc_type lu -richards_pc_factor_mat_solver_type superlu_dist -tadr_ksp_type preonly -tadr_pc_type lu -tadr_pc_factor_mat_solver_type superlu_dist"
```

`MMS_N` = 10, 20, 40, 80, 160, 320; `MMS_FCT=0` for low order; `MMS_DT` for dt/2
(the values are in `runs/jobs_dt2*.txt`).

## 3. Henry (`henry/` and `henry/low_order/`)

```bash
# henry/ (FCT) and henry/low_order/ (low order), the same command in each
RICHARDS_EDGE_VELOCITY=1 mpiexec -np 4 parun ls_CCS_so.py -l3 \
  -P "-ksp_type preonly -pc_type lu -pc_factor_mat_solver_type superlu_dist -tadr_ksp_type preonly -tadr_pc_type lu -tadr_pc_factor_mat_solver_type superlu_dist"
```

## 4. Elder (`elder/`)

```bash
RICHARDS_EDGE_VELOCITY=1 mpiexec -np 64 parun ls_CCS_so.py -l3 \
  -P "-richards_ksp_type preonly -richards_pc_type lu -richards_pc_factor_mat_solver_type superlu_dist -tadr_ksp_type preonly -tadr_pc_type lu -tadr_pc_factor_mat_solver_type superlu_dist"
```

## 5. Bioswale, bio2d set-up (`bioswale/`)

```bash
RICHARDS_EDGE_VELOCITY=1 mpiexec -np 2 parun ls_CCS_so.py -l3 \
  -P "-ksp_type preonly -pc_type lu -pc_factor_mat_solver_type superlu_dist -tadr_ksp_type preonly -tadr_pc_type lu -tadr_pc_factor_mat_solver_type superlu_dist"
```

## 6. Levee (`levee_density/density/` and `levee_density/non_density/`)

```bash
# levee_density/density/ (rho_s = 1.025) and levee_density/non_density/ (passive tracer)
RICHARDS_EDGE_VELOCITY=1 mpiexec -np 64 parun ls_CCS_so.py -l3 \
  -P "-ksp_type preonly -pc_type lu -pc_factor_mat_solver_type superlu_dist -tadr_ksp_type preonly -tadr_pc_type lu -tadr_pc_factor_mat_solver_type superlu_dist"
```
