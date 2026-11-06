#!/usr/bin/env bash
# usage: run.sh [run dir]   (default: ./run)
#
# Levee cross-section, x = 0..160 m, z = 0..57 m, homogeneous sand, density-coupled
# vs passive-tracer runs.  Decks from Levee_contamination_Density/Case_D and
# Case_D/non_density (commit 5066f6b5), with TADR switched to STAB 5.
#   flow    Richards STAB 2 + FCT, IC psi = 45 - z (water table 45 m),
#           psi = 48 - z on the wetted river face (leftTop), no flow on the
#           rest of the left side, top and bottom, and no flux condition set
#           on the right (natural outflow)
#   trans   TADR STAB 5 (ImplicitEV, backward Euler at CFL 0.5) + FCT, advected
#           with Richards' edge velocity; c_in = 0.7 on the wetted river face,
#           Dm = 1.886e-2 m^2/d, alpha_L = 0.1 m, alpha_T = 0.01 m
#   density density/: rho_s = 1.025 (eps = 0.025); non_density/: rho_s = 1,
#           DENSITY_MODEL=None.  Nothing else differs (a0 differs but is unused)
#   time    T = 2.0e4 d, nDTout = 801; the plots show archived steps 100, 400
#           and the last: t = 2446.9, 9937.6, 20025.0 d
# Both runs go at once, LEVEE_NP ranks each (64).  Each writes u_range.txt, the
# u range over the whole run, to check against the 0.7 inflow.  Solver options:
# petsc.options.superlu.  Detach, which otherwise dies with the session:
#   setsid nohup ./run.sh > run.log 2>&1 < /dev/null & disown
set -u
DECK="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUN="$(mkdir -p "${1:-$DECK/run}" && cd "${1:-$DECK/run}" && pwd)"

source ~/miniforge3/etc/profile.d/conda.sh
conda activate tadr-check
# parun must import the installed build, not a repo checkout on PYTHONPATH
unset PYTHONPATH
export RICHARDS_EDGE_VELOCITY=1
export OMP_NUM_THREADS=1

for v in density non_density; do
  (
  mkdir -p "$RUN/$v"
  cp "$DECK/$v"/*.py "$DECK"/petsc.options.superlu "$DECK/../vdft_metrics.py" "$RUN/$v"/
  cd "$RUN/$v"
  {
    echo "date: $(date)"
    echo "conda env: $CONDA_DEFAULT_ENV"
    python -c "import proteus,os;print('proteus:', os.path.dirname(proteus.__file__))"
    echo "repo HEAD: $(git -C "$DECK" rev-parse HEAD)"
    echo "np: ${LEVEE_NP:-64}"
  } > provenance.txt
  # uncommitted source edits at run time; the installed build may still differ
  git -C "$(git -C "$DECK" rev-parse --show-toplevel)" diff HEAD -- proteus > source.diff
  mpiexec -np "${LEVEE_NP:-64}" parun ls_CCS_so.py -l3 -O petsc.options.superlu > run.out 2>&1
  echo "exit $?" >> run.out
  # u range over every step, from the [Coupling rho q] lines; max must stay <= 0.7
  grep -o 'u (min,max,mean)=([^)]*' ls_CCS_so.log | tr '(,' '  ' | \
    awk '{if (NR==1 || $5<lo) lo=$5; if (NR==1 || $6>hi) hi=$6} END {print "u range over the run:", lo, hi}' > u_range.txt
  python plot_u_three_times_overlay.py && mv u_three_times_overlay.png "u_three_times_overlay_$v.png"
  ) &
done
wait

cd "$RUN"
python "$DECK/plot_domain_levee.py"    # levee_river_side_no_plume.png
