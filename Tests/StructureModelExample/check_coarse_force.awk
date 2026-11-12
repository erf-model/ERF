# Checks the "structure coarse source:" lines ERF prints with erf.v > 0: the force the source
# restricted onto a coarser level puts on the air must equal the force of the model's source
# on its own level, per direction, to tol (relative to the largest component). It fails when
# there is no such line or no force at all.
#   awk -v tol=1e-12 -f check_coarse_force.awk run.log
/^structure coarse source: level / {
    n++
    scale = 0
    for (d = 0; d < 3; d++) { v = $(15 + d); if (v < 0) v = -v; if (v > scale) scale = v }
    if (scale == 0) { zero++; next }
    for (d = 0; d < 3; d++) {
        diff = $(8 + d) - $(15 + d); if (diff < 0) diff = -diff
        if (diff > tol * scale) { bad++; print "force differs: " $0 }
    }
}
END {
    print n " coarse source lines, " bad + 0 " differ, " zero + 0 " without a force"
    exit !(n > 0 && bad == 0 && zero < n)
}
