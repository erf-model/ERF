# POSIX awk: use floating point for the area-weighted budget oracle (CMake's
# math(EXPR) is integer-only). Each argument is <leg>/cloud_chamber_budget.dat.
function absval(x) { return x < 0 ? -x : x }
function fail(msg) { print "CloudChamber_SidewallBands: " msg; failed = 1 }
function near_val(label, got, want, tol) {
    if (!(absval(got - want) <= tol))
        fail(sprintf("%s: got %.12g, expected %.12g +/- %.3g", label, got, want, tol))
}
BEGIN {
    split("cold hot one_band panels", legs)
    split("rhoTheta water_vapor", scalars)
    split("xlo xhi ylo yhi zlo zhi", faces)
    number = "^[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][+-]?[0-9]+)?$"
}
FNR == 1 {
    split(FILENAME, path, "/"); leg = path[1]
    if ($1 != "#" || $6 != "scalar" || $8 != "xlo" || $11 != "yhi" || $15 != "volume_change")
        fail(FILENAME ": unrecognized budget header")
    next
}
$5 == "rhoTheta" || $5 == "water_vapor" {
    s = $5; seen[leg,s]++
    if (NF != 18 || $1 != 0 || $2 != 0 || $3 != 1 || $4 !~ number || absval($4 - 0.0001) > 1e-10 ||
        $18 != (s == "water_vapor" ? "PASS" : "UNSUPPORTED_SOURCE"))
        fail(FILENAME ": unexpected step, time or status: " $0)
    for (i = 7; i <= 17; i++)
        if ($i !~ number) fail(FILENAME ": nonfinite or invalid budget field: " $i)
    for (f = 1; f <= 6; f++) {
        flux[leg,s,f] = $(6+f) + 0
        if (f != 1 && f != 4 && !(s == "rhoTheta" && f == 2))
            near_val(leg " " s " insulated/impermeable " faces[f], flux[leg,s,f], 0, 0)
    }
    if (s == "water_vapor") {
        # Coordinate flux is positive along the axis, so evaporating yhi is
        # negative. The state starts dry and remains far below saturation:
        # both incoming fluxes must actually reach the conserved vapor state.
        gain = $7 - $10
        near_val(leg " vapor volume_change", $14, gain, 5e-4 * absval(gain))
    }
}
END {
    for (l in legs) for (k in scalars)
        if (seen[legs[l],scalars[k]] != 1) fail(legs[l] ": need one " scalars[k] " row")
    if (failed) exit 1

    for (k in scalars) for (f = 1; f <= 6; f++) {
        s = scalars[k]; label = s " " faces[f]
        cold = flux["cold",s,f]; hot = flux["hot",s,f]; panels = flux["panels",s,f]
        near_val("one full band vs scalar " label, flux["one_band",s,f], cold, 2e-5 * absval(cold))
        if (f == 1 || f == 4) {
            sign_in = (f == 1 ? 1 : -1)
            if (s == "water_vapor" && !(sign_in*cold > 1e-9 && sign_in*hot > 2*sign_in*cold))
                fail(label ": scalar controls did not resolve cold and hot evaporation")
            if (s == "rhoTheta" && !(sign_in*cold < 0 && sign_in*hot > 0))
                fail(label ": scalar controls did not resolve cooling and heating")
            # Wall width and height are both 1 m. Hot areas are 0.90625-0.1875
            # on xlo and 0.265625 on yhi. Mix measured fluxes, not temperatures.
            hot_fraction = (f == 1 ? 23/32 : 17/64)
            # This is NOT exact superposition of RK2 solutions. Across 1 m,
            # rho*g*H/p ~ 1.2e-4; the first-stage wall exchange changes the air
            # by ~2*alpha*dt/dx^2 = 1.28e-5 of the contrast (twice at corners).
            # 5e-4 of the scalar span covers both, far below a lost panel or
            # evaluating saturation at the wall-mean temperature.
            near_val("panels " label, panels, cold + hot_fraction*(hot-cold), 5e-4 * absval(hot-cold))
        } else if (s == "rhoTheta" && f == 2) {
            if (!(cold < 0)) fail("scalar xhi control should heat the air")
            near_val("panels unchanged scalar xhi heat", panels, cold, 5e-4 * absval(cold))
        }
    }
    exit failed ? 1 : 0
}
