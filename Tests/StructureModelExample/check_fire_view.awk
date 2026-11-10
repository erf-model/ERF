# Checks the "structure example fire:" lines the example structure model prints when
# the run has the example fire engine (Tests/FireEngineExample): a disk of radius
# r(t) = r0 + ros t releasing q0 [W/m2] inside it.
#
# Limits: the engine must take one substep per step (fire.example.nsub = 1; with more,
# its flux is the mean over the substeps, not the flux at the step's midpoint); every
# step of the fire's level has the same length dtf, and the coarse steps start at
# multiples of lead (a fixed step from t = 0, also after a restart from such a time).
#
# Variables (awk -v):
#   lev     the level the fire grid refines
#   dtf     the step of that level [s]
#   lead    0 when the fire is on the model's level or a finer one: its fields belong to
#           the fire step that ended at the model's step time, and the two printed times
#           (17 digits) must agree exactly; otherwise the step of the
#           fire's (coarser) level [s]: the fire has already stepped to the end of the
#           coarse step that contains the model's step
#   r0, ros, q0  the example engine's inputs
#   nlines  the number of model steps expected
#   tol     relative tolerance of the heat and area against the disk (fire cells of 10 m
#           resolve the heat of these disks to 1.4 per cent or better, wherever they sit;
#           the example engine's burned_area() is exact)
#
# Expected, with tf the time of the fire's fields: tf = t, or the end of the coarse step;
# heat = q0 pi r(tf - dtf/2)^2 (the flux of the step that ended at tf; 0 before the
# fire's first step at t = 0, while after a restart the engine restores the flux of
# the step before it); burned = pi r(tf)^2.
function abs(x) { return x < 0 ? -x : x }
BEGIN { pi = atan2(0, -1); n = 0; bad = 0 }
/^structure example fire: / {
    n++
    t = $6 + 0; flev = $9 + 0; tf = $12 + 0; heat = $15 + 0; burned = $18 + 0
    want_tf = (lead > 0) ? lead * (int(t / lead + 1.0e-9) + 1) : t
    if (flev != lev)                   { print "line " n ": fire_lev = " flev ", expected " lev; bad = 1 }
    if (lead > 0) {
        if (abs(tf - want_tf) > 1.0e-9) { print "line " n ": fire_t = " tf " at t = " t ", expected " want_tf; bad = 1 }
    } else if ($12 != $6)              { print "line " n ": fire_t = " $12 ", expected the step time " $6 " exactly"; bad = 1 }
    r = r0 + ros * tf
    want_area = pi * r * r
    if (abs(burned - want_area) > tol * want_area) {
        print "line " n ": burned = " burned " at fire_t = " tf ", expected " want_area; bad = 1
    }
    if (tf == 0) {
        if (heat != 0) { print "line " n ": heat = " heat " before the fire's first step, expected 0"; bad = 1 }
    } else {
        rm = r0 + ros * (tf - 0.5 * dtf)
        want_heat = q0 * pi * rm * rm
        if (abs(heat - want_heat) > tol * want_heat) {
            print "line " n ": heat = " heat " at fire_t = " tf ", expected " want_heat; bad = 1
        }
    }
}
END {
    if (n != nlines) { print n " fire lines, expected " nlines; bad = 1 }
    if (!bad) { print "fire view checked on " n " steps" }
    exit bad
}
