#!/bin/sh
# The momentum budget and the restart of the reflux of the pressure impulse on coarse/fine
# faces (erf.momentum_reflux, erf-model/ERF#5017), on the example structure model: one point
# pushes on still air with F = 50 N along x from a refined level 1 (ratio 2, two-way coupling).
#
# usage: momentum_reflux_budget.sh "<mpi launcher>" <erf_exec> <fvolumesum> <fcompare>
#
# Run from a directory that holds inputs_structure_example. It checks, for a patch whose
# coarse/fine faces are normal to the force (x 4-11) and for a cube (x, y, z 4-11):
#   1. the composite integral of x_velocity at t = 4 s plus the impulse the reflux has still to
#      give (the last "Momentum reflux at level 0" line, x part, over rho_0 = 1) is
#      F t / rho_0 = 200 m^4/s to 1e-6 relative (it is 4e-9; with the option off 3e-3 and 1e-2);
#   2. a restart from the mid-run checkpoint reads that impulse and gives the same final
#      plotfile, bit for bit, and the same last impulse;
# that a checkpoint written with the option off restarts with it on (no impulse to read);
# and that a level whose fine level goes away gives its last impulse once and then none: with
# the model on level 0 and a level 1 away from it, tagged where u < 1e-12 m/s, which the
# regrid before the third step removes, the integral grows by F dt / rho_0 = 12.5 m^4/s per
# step from the end of the third step on (162.5 m^4/s from t = 0.75 s to t = 4 s, to 1e-6
# relative). An impulse given in the fourth step too, or in every step, fails this.
set -e
MPI=$1
EXE=$2
FVS=$3
FCMP=$4

COMMON="erf.plot_int_1=16 structure.example.epsilon=25 amr.max_level=1 erf.dt_ref_ratio=2"
COMMON="$COMMON erf.coupling_type=TwoWay structure.example.level=1 erf.refinement_indicators=box1"
COMMON="$COMMON erf.box1.max_level=1 erf.v=1 amrex.call_addr2line=0"

run () {  # the launcher is quoted for the shell (it may hold a path with spaces)
    eval "$MPI \"\$EXE\" $*"
}

pending () {  # log: the x part of the last impulse still to be given on level 0
    awk '/Momentum reflux at level 0:/ { for (i = 1; i < NF; i++) { if ($i == "x") { v = $(i+2) + 0 } } }
         END { if (v == "") { exit 1 } printf "%.16e\n", v }' "$1"
}

one_case () {  # name, low and high corner of the patch in level-0 cells
    name=$1
    box="erf.box1.in_box_lo_indices_crse=$2 erf.box1.in_box_hi_indices_crse=$3"
    rm -rf "$name"
    mkdir -p "$name/rst"
    cd "$name"
    run ../inputs_structure_example $COMMON $box erf.momentum_reflux=true erf.check_int=8 > run.log 2>&1
    integral=$("$FVS" -v x_velocity plt00016 | awk '/integral of x_velocity/ { print $NF }')
    impulse=$(pending run.log)
    awk -v n="$name" -v b="$integral" -v p="$impulse" 'BEGIN {
        e = 200.0; d = b + p - e; if (d < 0) { d = -d }
        printf "%s: integral of x_velocity = %s, impulse still to be given = %s, sum - 200 = %.3e\n", n, b, p, b + p - e
        exit !(d < 1.0e-6 * e) }'
    cd rst
    run ../../inputs_structure_example $COMMON $box erf.momentum_reflux=true erf.check_int=-1 erf.restart=../chk00008 > ../restart.log 2>&1
    cd ..
    grep -q "Momentum reflux: read the impulse still to be given on level 0" restart.log
    "$FCMP" -r 0 --abs_tol 0 plt00016 rst/plt00016
    test "$(pending restart.log)" = "$impulse"
    cd ..
}

one_case normal "4 0 0" "11 15 15"
one_case cube   "4 4 4" "11 11 11"

# A checkpoint without the impulse (the option was off): the restart runs and reads none
rm -rf nofile
mkdir -p nofile/rst
cd nofile
run ../inputs_structure_example $COMMON erf.box1.in_box_lo_indices_crse=4 0 0 erf.box1.in_box_hi_indices_crse=11 15 15 max_step=8 erf.plot_int_1=-1 erf.check_int=8 > run.log 2>&1
test ! -e chk00008/Level_0/MomRefluxX_H
cd rst
run ../../inputs_structure_example $COMMON erf.box1.in_box_lo_indices_crse=4 0 0 erf.box1.in_box_hi_indices_crse=11 15 15 erf.momentum_reflux=true erf.plot_int_1=-1 erf.check_int=-1 erf.restart=../chk00008 > ../restart.log 2>&1
cd ..
if grep -q "Momentum reflux: read the impulse" restart.log; then exit 1; fi
grep -q "Momentum reflux at level 0:" restart.log
cd ..

# A level 1 that goes away: the impulse of the last reflux is given once, in the next step
rm -rf gone
mkdir gone
cd gone
run ../inputs_structure_example structure.example.epsilon=25 structure.example.level=0 amr.max_level=1 erf.dt_ref_ratio=2 erf.coupling_type=TwoWay erf.momentum_reflux=true erf.v=1 amrex.call_addr2line=0 erf.check_int=-1 erf.plot_int_1=1 erf.regrid_int=2 erf.refinement_indicators=r1 erf.r1.max_level=1 erf.r1.field_name=u erf.r1.value_less=1.e-12 erf.r1.in_box_lo=300 150 150 erf.r1.in_box_hi=400 250 250 > run.log 2>&1
test "$(grep -c 'Momentum reflux at level 0:' run.log)" = 2
# The check below sees an impulse given more than once only if that impulse is not small
awk -v p="$(pending run.log)" 'BEGIN { if (p < 0) { p = -p }; exit !(p > 1.0e-2) }'
i03=$("$FVS" -v x_velocity plt00003 | awk '/integral of x_velocity/ { print $NF }')
i16=$("$FVS" -v x_velocity plt00016 | awk '/integral of x_velocity/ { print $NF }')
awk -v a="$i03" -v b="$i16" 'BEGIN {
    e = 162.5; d = b - a - e; if (d < 0) { d = -d }
    printf "gone: integral of x_velocity at t = 0.75 s = %s, at t = 4 s = %s, difference - 162.5 = %.3e\n", a, b, b - a - e
    exit !(d < 1.0e-6 * e) }'
cd ..
