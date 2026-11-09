# Checks the "structure example fire:" lines the example structure model prints when
# the run has the example fire engine (Tests/FireEngineExample) with nsub = 1: a disk
# of radius r(t) = r0 + ros t releasing q0 [W/m2] inside it.
#
# Variables (awk -v):
#   lev     the level the fire grid refines
#   dtf     the step of that level [s]
#   lead    0 when the fire is on the model's level or a finer one: its fields belong to
#           the fire step that ended at the model's step time; otherwise the step of the
#           fire's (coarser) level [s]: the fire has already stepped to the end of the
#           coarse step that contains the model's step
#   r0, ros, q0  the example engine's inputs
#   nlines  the number of model steps expected
#   tol     relative tolerance of the heat and area against the disk (the fire cells
#           resolve the circle to about 1 per cent)
#
# Expected, with tf the time of the fire's fields: tf = t, or the end of the coarse step;
# heat = q0 pi r(tf - dtf/2)^2 (the flux of the step that ended at tf; 0 before the
# fire's first step at t = 0); burned = pi r(tf)^2.
function abs(x) { return x < 0 ? -x : x }
BEGIN { pi = atan2(0, -1); n = 0; bad = 0 }
/^structure example fire: / {
    n++
    t = $6 + 0; flev = $9 + 0; tf = $12 + 0; heat = $15 + 0; burned = $18 + 0
    want_tf = (lead > 0) ? lead * (int(t / lead + 1.0e-9) + 1) : t
    if (flev != lev)                   { print "line " n ": fire_lev = " flev ", expected " lev; bad = 1 }
    if (abs(tf - want_tf) > 1.0e-9)    { print "line " n ": fire_t = " tf " at t = " t ", expected " want_tf; bad = 1 }
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
