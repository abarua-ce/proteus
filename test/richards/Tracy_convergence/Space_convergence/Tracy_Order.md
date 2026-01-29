
##  Proteus branch 
```
git checkout abarua-ce/Tracy_convergence_ceds
```
```
git branch 
Tracy_convergence_ceds
cd ~/proteus/test/richards/Tracy_convergence/Space_convergence
ls
FCT  Stab_0  Stab_2
```


Each stabilization method contains multiple mesh refinement levels:

```
ref_0  ref_1  ref_2  ref_3  ref_4  ref_5
```

---

# Running the Simulations

## Step 1 — Enter a Stabilization Folder

Example:

```bash
cd FCT
```

Repeat the same workflow for:

```
Stab_0
Stab_2
```

---

## Step 2 — Run Mesh Refinements (Coarse → Fine)

Run the following commands inside each refinement directory.

---

### ref_0

```bash
cd ref_0

mpiexec -n 4 parun re_vgm_sand_10x10m_2d_p.py re_vgm_sand_10x10m_2d_c0p1_n.py -l 5-v  -P "-ksp_type preonly -pc_type lu -pc_factor_mat_solver_type superlu_dist"
```
---

### ref_1

```bash
cd ../ref_1

mpiexec -n 8 parun re_vgm_sand_10x10m_2d_p.py \
re_vgm_sand_10x10m_2d_c0p1_n.py \
-l 5 -v \
-P "-ksp_type preonly -pc_type lu -pc_factor_mat_solver_type superlu_dist"
```

---

### Continue for Higher Refinements

| Refinement | nnx | Suggested MPI Ranks |
|------------|------|--------------------|
| ref_0 | 11  | 4 |
| ref_1 | 21  | 8 |
| ref_2 | 41  | 8–16 |
| ref_3 | 81  | 16 |
| ref_4 | 161 | 32 |
| ref_5 | 321 | 32–64 |


---

# Create seperate environment for Error Computation that have quadpy-legacy

I have created seperate environment where I only calculate Convergence Order. The convergence analysis uses **quadpy-legacy** for high-accuracy numerical integration.

---

##  Create Conda Environment

Recommended Python version: **3.11**

```bash
conda create -n quadpy_env python=3.11 -y
conda activate quadpy_env
conda install scipy h5py -c conda-forge
python -m pip install legacy-quadpy==0.16.10
```
---



---

# Compute Order of Accuracy

After **all refinement simulations are finished**, navigate to the stabilization directory (example: `FCT`) and run:

```bash
python Error_rate.py
```

This script will:

- Read simulation outputs  
- Compute L2 and L∞ error norms  
- Estimate the spatial order of accuracy  
- Generate convergence tables  

---

---

# Quick Workflow Summary

```bash
git checkout Tracy_convergence_ceds

cd ~/proteus/test/richards/Tracy_convergence/Space_convergence/FCT

# Run refinements
cd ref_0 -> run
cd ../ref_1 -> run
...
cd ../ref_5 -> run

# Compute convergence
conda activate quadpy_env
python Error_rate.py
Do same for Stab_2 and Stab_0  
```

---
