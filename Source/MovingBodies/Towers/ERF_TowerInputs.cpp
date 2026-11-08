// TowerType::force_coefficient() and the range checks of a tower type's inputs.

#include "ERF_TowerInputs.H"

#include <cmath>
#include <utility>

#include "ERF_MemberDrag.H"

using amrex::Real;

namespace erf_towers {

namespace {
// finite and > 0; finite and >= 0 (NaN and infinity fail both)
bool positive (Real x) { return std::isfinite(x) && x > 0.0; }
bool non_negative (Real x) { return std::isfinite(x) && x >= 0.0; }
}

Real TowerType::force_coefficient () const
{
    return drag_coefficient > 0.0 ? drag_coefficient : lattice_force_coefficient(solidity);
}

std::string TowerType::validate () const
{
    const std::string key = "erf.conductors." + name + ".";
    if (!positive(base_width)) { return key + "base_width must be finite and positive (m)"; }
    if (!(positive(top_width) && top_width <= base_width)) { return key + "top_width must be finite, positive and at most base_width (m)"; }
    if (!(solidity > 0.0 && solidity < 1.0)) { return key + "solidity must be in (0, 1): the members' area over a face's outline"; }
    if (!positive(arm_length)) { return key + "arm_length must be finite and positive (m)"; }
    if (arm_outside_shaft && !(arm_length > top_width)) {
        // the least it can be: each tower also needs it longer than the shaft is wide at the middle of the arm's face
        return key + "arm_length must exceed top_width for " + key + "arm_outside_shaft (the cross-arm's drag outside the "
                     "shaft); set " + key + "arm_outside_shaft = false for a cross-arm within the shaft's width";
    }
    if (!non_negative(arm_depth)) { return key + "arm_depth must be finite and >= 0 (m; 0: top_width)"; }
    if (!non_negative(peak)) { return key + "peak must be finite and >= 0 (m)"; }
    if (!non_negative(drag_coefficient)) { return key + "drag_coefficient must be finite and >= 0 (0: from the solidity)"; }
    if (segments < 1) { return key + "segments must be >= 1"; }
    if (!std::isfinite(epsilon)) { return key + "epsilon must be finite"; }
    if (epsilon_given && !(epsilon > 0.0)) {
        return key + "epsilon must be positive (in units of the anchor level's dx; not given: erf.conductors.epsilon)";
    }
    if (!non_negative(weight)) { return key + "weight must be finite and >= 0 (N)"; }
    if (!non_negative(leg_spacing)) { return key + "leg_spacing must be finite and >= 0 (m; 0: base_width)"; }
    if (!non_negative(allowable_uplift)) { return key + "allowable_uplift must be finite and >= 0 (N; 0: not checked)"; }
    if (!non_negative(allowable_compression)) { return key + "allowable_compression must be finite and >= 0 (N; 0: not checked)"; }
    if (!non_negative(frequency)) { return key + "frequency must be finite and >= 0 (Hz; 0: the tower stands still)"; }
    if (frame_panels < 0 || frame_panels > 200) { return key + "frame_panels must be in [0, 200] (0: no generated frame)"; }
    if (!frame_file.empty() && frame_panels > 0) {
        return key + "frame_file and " + key + "frame_panels both give the tower's frame; give one";
    }
    if (has_frame()) {
        // the frame model gives the stiffness, the mass and the footings
        const std::string source = key + (frame_file.empty() ? "frame_panels" : "frame_file");
        const std::pair<const char*, Real> from_frame[] = {
            {"frequency", frequency}, {"weight", weight}, {"foundation_rotational_stiffness", foundation_rotational_stiffness},
            {"foundation_lateral_stiffness", foundation_lateral_stiffness}};
        for (const auto& kv : from_frame) {
            if (kv.second != 0.0) {
                return key + kv.first + " is not given with " + source + ": the frame model sets the tower's stiffness, "
                       "mass and footings";
            }
        }
    }
    if (!member_file.empty() && frame_file.empty()) {
        return key + "member_file needs " + key + "frame_file, whose members it describes (a generated frame has its own)";
    }
    const std::pair<const char*, const std::vector<Real>*> angles[] = {{"leg_angle", &leg_angle}, {"brace_angle", &brace_angle}};
    for (const auto& kv : angles) {
        const std::vector<Real>& a = *kv.second;
        if (frame_panels == 0) {
            if (!a.empty()) { return key + kv.first + " needs " + key + "frame_panels, the frame it sizes"; }
            continue;
        }
        if (a.size() != 2 || !positive(a[0]) || !positive(a[1]) || !(a[1] < a[0])) {
            return key + kv.first + " needs two values, the equal-leg angle's leg width and thickness (m), with 0 < thickness < width";
        }
    }
    if (frame_panels == 0 && !bracing.empty()) { return key + "bracing needs " + key + "frame_panels, the frame it braces"; }
    if (!bracing.empty() && bracing != "crossed" && bracing != "single") {
        return key + "bracing must be crossed or single, not '" + bracing + "'";
    }
    if (!non_negative(yield_strength)) { return key + "yield_strength must be finite and >= 0 (Pa; 0: 3.45e8)"; }
    if (frame_panels == 0 && yield_strength != 0.0) {
        return key + "yield_strength needs " + key + "frame_panels (the members of frame_file take theirs from member_file)";
    }
    if (!(std::isfinite(steel_temperature) && steel_temperature > -273.15 && steel_temperature < 1200.0)) {
        return key + "steel_temperature must be above -273.15 and below 1200 (C): the steel keeps no stiffness at 1200 C";
    }
    if (steel_temperature != 20.0 && !has_frame()) {
        return key + "steel_temperature needs a frame model (" + key + "frame_file or " + key + "frame_panels), whose steel it heats";
    }
    if (frequency > 0.0 && !(weight > 0.0)) { return key + "frequency needs the tower's weight (N), which sets its mass"; }
    if (!(damping_ratio >= 0.0 && damping_ratio < 1.0)) { return key + "damping_ratio must be in [0, 1) (fraction of critical)"; }
    if (!non_negative(foundation_rotational_stiffness)) {
        return key + "foundation_rotational_stiffness must be finite and >= 0 (N m/rad; 0: rigid)";
    }
    if (!non_negative(foundation_lateral_stiffness)) { return key + "foundation_lateral_stiffness must be finite and >= 0 (N/m; 0: rigid)"; }
    if (!moves() && !frequency_given) {
        // a tower that stands still has no motion for these to act on (with frequency = 0 given, a warning:
        // unused_on_a_still_tower)
        if (damping_given) { return key + "damping_ratio needs a tower that moves: " + key + "frequency, frame_file or frame_panels"; }
        const std::pair<const char*, bool> footing[] = {
            {"foundation_rotational_stiffness", foundation_rotational_stiffness != 0.0},
            {"foundation_lateral_stiffness", foundation_lateral_stiffness != 0.0}};
        for (const auto& kv : footing) {
            // a frame's supports are its footings: only a one-mode tower takes these
            if (kv.second) { return key + kv.first + " needs a tower that moves in one mode: " + key + "frequency"; }
        }
    }
    if (has_frame() && !(damping_ratio > 0.0)) {
        // the members' drag is held from the start of each coupling step and handed to Newmark's method at its
        // end: a damping force one step late, which an undamped frame's higher modes (from 6.4 times the first
        // natural frequency at 20 coupling steps a period) grow under
        return key + "damping_ratio must be positive for a frame (frame_file or frame_panels): without structural "
                     "damping the members' drag makes its higher modes grow";
    }
    if (has_frame() && leg_spacing != 0.0) {
        return key + "leg_spacing is not given with a frame: the frame's supports are its footings";
    }
    if (angle_axes_given && frame_panels == 0) {
        return key + "angle_principal_axes needs " + key + "frame_panels (a generated frame; a frame_file sets each member's MSpin)";
    }
    return std::string();
}

std::string TowerType::unused_on_a_still_tower () const
{
    if (moves() || !frequency_given) { return std::string(); }
    std::string keys;
    auto add = [&] (const char* k) { keys += (keys.empty() ? "" : ", ") + std::string(k); };
    if (damping_given) { add("damping_ratio"); }
    if (foundation_rotational_stiffness != 0.0) { add("foundation_rotational_stiffness"); }
    if (foundation_lateral_stiffness != 0.0) { add("foundation_lateral_stiffness"); }
    if (keys.empty()) { return keys; }
    return "erf.conductors." + name + ": " + keys + " not used: frequency = 0 and no frame, so the tower stands still";
}

} // namespace erf_towers
