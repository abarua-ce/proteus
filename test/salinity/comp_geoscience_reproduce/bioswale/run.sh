#!/usr/bin/env bash
# usage: run.sh [run dir]   (default: ./run)
#
# Bioswale cross-section, 3 m x 5 m, set up as test/richards/WRR_input_files/bio2d
# except for the drain: here a 0.25 m-radius circle (domain_rg.py, nnx = 21),
# there a rectangle.  Like bio2d the drain wall is no-flux (useSeepageFace = False).
#   flow    Richards STAB 2 + FCT, uniform psi = -3 m IC, 0.5 m ponding on the top
#           for 1 <= x <= 2 m, three van Genuchten zones; Newton/step control as
#           in bio2d (floor 10, ceil 20, grow 2, 5 line searches)
#   trans   TADR STAB 5 (ImplicitEV, backward Euler at CFL 0.5) + FCT,
#           c_in = 1, Dm = 1.886e-3, alpha_L = 0.1, alpha_T = 0.01, rho_s = 1.025
#   time    T = 0.5 d, nDTout = 201; the plot shows t = 0.1219, 0.2488, 0.5025 d
# 2 ranks by default (BIOSWALE_NP).  Solver options: petsc.options.superlu.
set -u
DECK="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUN="$(mkdir -p "${1:-$DECK/run}" && cd "${1:-$DECK/run}" && pwd)"
cp "$DECK"/{domain_rg,Raingarden_p,Raingarden_n,tranport_p,tranport_n,ls_CCS_so,plot_u_three_times,plot_u_velocity,plot_domain_bioswale}.py "$DECK"/petsc.options.superlu "$DECK/../vdft_metrics.py" "$RUN"/
cd "$RUN"

source ~/miniforge3/etc/profile.d/conda.sh
conda activate tadr-check
# parun must import the installed build, not a repo checkout on PYTHONPATH
unset PYTHONPATH
export RICHARDS_EDGE_VELOCITY=1
export OMP_NUM_THREADS=1

{
  echo "date: $(date)"
  echo "conda env: $CONDA_DEFAULT_ENV"
  python -c "import proteus,os;print('proteus:', os.path.dirname(proteus.__file__))"
  echo "repo HEAD: $(git -C "$DECK" rev-parse HEAD)"
  echo "np: ${BIOSWALE_NP:-2}"
  echo "BIOSWALE_DENSITY: ${BIOSWALE_DENSITY:-0}"
} > provenance.txt
# uncommitted source edits at run time; the installed build may still differ
git -C "$(git -C "$DECK" rev-parse --show-toplevel)" diff HEAD -- proteus > source.diff

mpiexec -np "${BIOSWALE_NP:-2}" parun ls_CCS_so.py -l3 -O petsc.options.superlu > run.out 2>&1
echo "exit $?" >> run.out

python plot_u_three_times.py && mv u_three_times.png u_three_times_bioswale.png
python plot_u_velocity.py         # u_colour_velocity_bioswale.png: u, isochlors, Darcy streamlines
python plot_domain_bioswale.py    # bioswale_contamination_domain.png
