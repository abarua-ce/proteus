# Reproducing the figures

One folder per manuscript figure.  Each holds the exact deck that produced the
figure (copied from the run directory that wrote the PNG), the solver options,
the input data, the plot script, the figure itself, and a `run.sh` that reruns
the deck and redraws the figure.

| Folder | Figure(s) | Copied from | `run.sh` |
|---|---|---|---|
| `zalesak/` | `disc_r.png` (not kept; see below) | `test/TADR/thelper_tadr*.py` | serial, 2 runs |
| `mms/` | `paper_figures/` (tables `tab:mms_t3_*`, `mms_*.png`) | `test/salinity/mms_coupled` | 67 serial runs, 10 at a time |
| `henry/` | `iso_combined_FCT_ref_1.png`, `low_order/iso_combined_Stab_2_lo.png`, `domain_henry_minimal.png` | `Henry_problem/Henry_corrected/right_dbc/ref_1` | 64 ranks, 2 runs |
| `elder/` | `elder_u_colour_velocity_20yr.png`, `elder_isochlors_velocity_20yr.png` | `Elder_problem/stab_5/run_uhdcfg_20yr_nocap` | 64 ranks, ~10 min |
| `simmons_uhd/` | `fig8_uhd_colormap.png` | `benchmark_problem/UMD/run_uhd_l4_ramp` | 90 ranks, ~70 min |

Every `run.sh` takes an optional run directory (default `./run`, `./run_uhd`
for Simmons) and copies the deck into it, so the bundle itself stays clean.
Detach the long ones:

    setsid nohup ./run.sh > run.log 2>&1 < /dev/null & disown

## Build and provenance

- Every run uses the `tadr-check` conda env with `PYTHONPATH` unset, so `parun`
  imports proteus from that env's site-packages, not from this checkout.
  Rebuild and install before running if the source has changed.
- Each run writes `provenance.txt` (date, env, the proteus path actually
  imported, repo HEAD, the `SIMMONS_*`/`RICHARDS_*`/`TADR_*` knobs) and
  `source.diff` (uncommitted edits under `proteus/` at run time).
- The figures here were NOT made from a commit.  All of them ran on builds of
  the uncommitted working tree of branch `WRR_ifem_merge` (HEAD `012b97d1`
  when this folder was made).  The Elder (2026-10-29) and UHD (2026-10-27)
  runs predate the 2026-11-02 source changes (Richards FCT `d^G` term removed;
  TADR FCT bounds on `c^{n+1}`), so today's build will not reproduce them bit
  for bit.  Henry ref_1 ran in May 2026 in the `petsc-dev-again` env.
- Solver settings: each deck's `petsc.options.superlu` (direct SuperLU_dist LU
  for both `richards_` and `tadr_`).  Henry passes the same options
  (`-ksp_type preonly -pc_type lu -pc_factor_mat_solver_type superlu_dist`) on
  the command line, as its original run did.

## Per case

### Zalesak's disk (`zalesak/`)
`thelper_tadr.py` with `problem=1`: notched disk, cone and hump on the unit
square, rigid rotation `v = 2 pi (0.5 - y, x - 0.5)`, so `T = 1` is one
revolution.  `refinement=4`: 161 x 161 nodes; STABILIZATION_TYPE 2
(EntropyViscosity), explicit SSP2, `cfl = 0.5`.  Low order is the same scheme
with `FCT=False` (TADR then advances with `uLow`).  `plot_disc.py` draws exact
/ low order / FCT and prints `u_min`, `u_max` and the L1 error of each.
**These settings were chosen on 2026-11-03.  The run that made the manuscript's
`disc_r.png` was not recorded, and `plot_disc.py` is new.  Rerun and replace
the figure.**

### Manufactured solution (`mms/`)
Definitions of T1-T3 and the two refinement paths: `MMS_PLAN.md` and
`mms_coupled_subsection.tex`.  Exact fields `psi_e`, `c_e` and the sources
`f_w` (flow) and `g` (transport): `mms_exact.py`.  Van Genuchten parameters
for the missing `tab:vgm` (UHD sand, `mms_exact.py`):
`K0 = 1.12e-3 m/s`, `alpha = 11.6 1/m`, `n = 2.39`, `m = 1 - 1/n`,
`theta_r = 0.0162`, `theta_s - theta_r = 0.3438`, plus `Dm = 1e-6`,
`eps = 0.235`.  `run.sh` runs `runs/jobs_*.txt` through `run_mms.sh` and writes
the tables and figures into `figures/`; `paper_figures/` keeps the 2026-11-02
results.
**Not covered: there are no T0 runs (Richards alone, `MMS_FLOW_ONLY=1`) and no
path (b) runs (steady, `MMS_STEADY=1`).  Path (a) was run with `dt` and `dt/2`
instead.**

### Henry (`henry/`)
The manuscript figure is **ref_1, not ref_3**: 2 m x 1 m, `nnx x nny = 80 x 40`
nodes (79 x 39 cells), `T = 3.0e4 s`, 201 outputs, TADR STAB 2, FCT on
(`low_order/` differs only in `FCT=False`).  ref_3 (320 x 160 nodes,
`T = 5e4 s`, STAB 4) is a later deck that no figure uses.
Analytical comparison: `henry_semianalytical.py` is the Henry semi-analytical
solution of Zidane et al. (2012, WRR 48, W06533, Appendix A), Henry case
`a = 0.263`, `b = 0.1`, `xi = 2`, Fourier truncation MA/NA/MB/NB = 7/4/4/19.
`run.sh` evaluates it at the run's nodes, then draws the 0.1/0.5/0.9 isochlors
of each scheme against it (`plot_iso_combined.py`), plus the full-field c
comparison with its RMS and max errors (`plot_henry_comparison.py`,
`henry_error_{fct,lo}.txt`).  `henry_reference.npz` is the reference the
figures used.

### Elder (`elder/`)
601 x 151 nodes, structured triangles, 600 m x 150 m; source on the top centre
(`tranport_p.py`, `conc_x0 = 150`); `p = 0` pinned at the two upper corners
(`flow_p.py`).
Saturated VG with `alpha = 1e-5`; `K = 0.4107 m/d`; constant viscosity.
Richards and TADR FCT on, TADR advected with Richards' edge velocity.
**Run length is 20 years** (`T = 7305 d`), output every 30.4375 d (240
outputs); the figures show 2, 4, 10 and 20 yr (Voss & Souza 1987 Fig. 9a).
The manuscript text's "4 years" should read 20 years.
`flow_n_base.py` and `tranport_n_base.py` are `Elder_problem/stab_5/flow_n.py`
and `tranport_n.py`, which the original run directory loaded as `../*.py`.

### Simmons UHD (`simmons_uhd/`)
`SIMMONS_LEVEL=4`: 500 x 500 elements.  Seepage drains at the two top corners,
time-varying inflow `data/q_UHD_ramp.csv` (digitised record `data/q_UHD.csv`
with a linear ramp over the unrecorded first 3 min), variable viscosity
(`SIMMONS_VISC=1`, `omega_max = 0.2535`), TADR STAB 5 + FCT.  Random ln K:
`k0_field.py`, sigma 0.01, seed 0, one value per grid rectangle; the
realisation is `data/k0_n250000_sig0.01_seed0_ln.npy`.  Guevara Morel & Graf
(2023) reference data: `data/z_hgs_UHD.csv` (HGS 30-realisation mean +/-
sigma) and `data/z_lab_UHD.csv` (lab), from `digitize_fig7.py`; lab-photo
masks at 32 and 40 min from `digitize_uhd_photos.py`.
**`fig8_uhd_colormap.png` shows t = 32 and 40 min** (`plot_fig8.py`,
`FIG_TIMES`).  The caption's 12/32/48 min is wrong.
