# Computational Geosciences runs

Decks and scripts for the figures and tables of the VDFT manuscript. Richards
(flow) STABILIZATION_TYPE 2 + FCT, TADR (transport) STABILIZATION_TYPE 5 + FCT,
coupled through the flow model's edge fluxes (`RICHARDS_EDGE_VELOCITY=1`, set by
every script). The low-order runs turn FCT off in both models.

## Conventions

- `<case>/run.sh [run dir]` copies the deck into the run directory (default
  `./run`), runs it and makes the plots. Settings are environment variables,
  listed in the header of each script.
- Long runs must be detached, or they die with the session:
  `setsid nohup ./run.sh > run.log 2>&1 < /dev/null & disown`
- Runs that copy `vdft_metrics.py` (bioswale, Elder, UHD) write
  `vdft_metrics.json`: steps, Newton iterations of both models, wall time, the
  solute mass-balance error, the drift of the carried fluid mass, the flow balance
  and the c range.
- Process counts below are those of `tab:solver`.

## 1. Zalesak's disk (`zalesak/`)

```bash
ZALESAK_REF=5 ZALESAK_NP=32 ./run.sh      # <run dir>/fct and <run dir>/low
```

Refinement 5 is the paper mesh (321 x 321 = 103,041 nodes, h = 1/320); the
script's defaults are refinement 4 and one process. `ZALESAK_STAB=2` runs the
explicit scheme instead.

## 2. Manufactured solution (`mms/`)

```bash
./run_mms.sh runs/a_T3_fct1_N40 MMS_N=40 MMS_FCT=1      # one run
```

- `MMS_N` = 10, 20, 40, 80, 160, 320; `MMS_FCT=0` for low order; `MMS_NP` ranks (1).
- Time-step studies: `MMS_DT` values in `runs/jobs_dt2.txt` (dt/2) and
  `runs/jobs_dth2.txt` (dt proportional to h^2, N = 10 to 80); each line is
  `<run dir> VAR=value ...`, the arguments of `run_mms.sh`.
- `run.sh` runs the convergence ladder; `mms_table.py` and `cost_table.py` build
  the tables.

## 3. Henry (`henry/`)

```bash
HENRY_NP=4 ./run.sh       # <run dir>/fct and <run dir>/lo (low_order/)
```

The errors against the semi-analytical solution come from `henry_metrics.py`,
which is on cekees/vdft-manuscript and not yet in this directory.

## 4. Elder (`elder/`)

```bash
ELDER_NP=64 ./run.sh run_dK                                   # paper run, adaptive dt
ELDER_DT=7.609375  ELDER_NP=64 ./run.sh run_dt7p6             # fixed dt = DT_out/4
ELDER_DT=3.8046875 ELDER_NP=64 ./run.sh run_dt3p8             # fixed dt = DT_out/8
ELDER_NNX=301 ELDER_NNY=76 ELDER_NP=64 ./run.sh run_coarse    # 301 x 76 mesh
```

`time_refinement_table.py` and `mesh_comparison_table.py` read these run
directories by name; `compare_dt.py` plots them.

## 5. Simmons UHD (`simmons_uhd/`)

```bash
SIMMONS_NP=64 ./run.sh run_uhd       # level 4, 500 x 500 cells, 48 min
```

- The run directory name must start with `run_uhd` (`plot_fig8.py` reads the
  case from it).
- Inflow: `data/q_UHD_ramp.csv`, the digitized record `data/q_UHD.csv` with q
  rising linearly from 0 at t = 0 to its first sample at 3 min. Reference depths:
  `data/z_lab_UHD.csv` (laboratory) and `data/z_hgs_UHD.csv` (HGS); the
  conductivity field is cached in `data/k0_*.npy` and regenerated from the seed
  if absent.
- Flow Newton tolerance `nl_atol_res = 1e-8` (`flow_n.py`); at 1e-7 the flow
  balance is about 1e-4 instead of 5e-6.

### Coupling ablation, level 0 (`tab:uhd_coupling`)

`run_level0.sh [run dir]` runs level 0 (10,201 nodes, 432 steps at Co 0.5) in
serial; the coupling is chosen by the environment:

| Environment | Coupling |
|---|---|
| unset | compatible: reconstructed flux, carried fluid mass |
| `TADR_VELOCITY_COUPLING=eos` | Galerkin Darcy velocity, storage theta(psi^n) rho(c^n) |
| `TADR_VELOCITY_COUPLING=lagged` | Galerkin Darcy velocity, storage of the preceding flow solve |
| `TADR_CONSERVATIVE_FORM=1` (with one of the above, `SIMMONS_FCT=0`) | conservative instead of cancelled solute row |

Also `SIMMONS_COURANT` (0.5), `SIMMONS_FCT` (1, transport limiter) and
`SIMMONS_RICHARDS_FCT` (1, flow limiter).

```bash
level0_rerun/run_all.sh          # the five rows: Co 0.5, Co 0.25, lagged, lagged with flow limiter off, compatible
level0_conservative/run_all.sh   # cancelled vs conservative form, transport limiter off
```

Each run directory's `simmons_so.log` holds the solute created by the coupling
(last `[TADR-CONT2]` line, flux and storage parts) and the final solute mass
(last `[TADR-MASS]` line).

## 6. Bioswale (`bioswale/`)

```bash
BIOSWALE_DENSITY=1 ./run.sh run_final     # density-coupled, drain wall impermeable
```

`BIOSWALE_DENSITY=0` runs the passive tracer, `BIOSWALE_SEEPAGE=1` makes the drain
wall a seepage face, `BIOSWALE_NNX` sets the mesh (21), `BIOSWALE_NP` the ranks (2).

## 7. Levee (`levee_density/`)

```bash
LEVEE_NP=4 ./run.sh      # density/ (rho_s = 1.025) and non_density/ (passive tracer), side by side
```

## Tables

`tables/make_tables.py` writes the Markdown and LaTeX tables in `tables/` from the
run directories it names; `bounds_mass_table.py LABEL:RUN_DIR:LO:HI ...` prints
the boundedness and mass rows for any set of runs.
