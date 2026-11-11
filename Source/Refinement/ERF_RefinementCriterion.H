#ifndef ERF_REFINEMENT_CRITERION_H_
#define ERF_REFINEMENT_CRITERION_H_

#include <AMReX_REAL.H>
#include <AMReX_RealBox.H>
#include <AMReX_Vector.H>

#include <algorithm>
#include <limits>
#include <string>

/**
 * One test in a refinement criterion: a field compared against a threshold that may differ by level
 */
struct TagCondition
{
    enum struct Test { Greater, Less, AdjacentDifferenceGreater };

    std::string field;
    Test test = Test::Greater;

    // Threshold at level lev is values[lev]; the last entry holds for all higher levels
    amrex::Vector<amrex::Real> values;

    [[nodiscard]] amrex::Real threshold (int lev) const {
        return values[std::min(lev, static_cast<int>(values.size())-1)];
    }
};

/**
 * A refinement indicator from erf.refinement_indicators, tagged by ERF itself.
 *
 * A cell is tagged when the conditions hold (all of them, or any of them, as require_all says)
 * and the cell is inside the box, if there is one. A criterion without conditions tags every
 * cell in its box.
 */
struct RefinementCriterion
{
    // Name in erf.refinement_indicators, which is also the ParmParse prefix erf.<name>
    std::string name;

    amrex::Vector<TagCondition> conditions;

    // Whether every condition must hold (erf.<name>.combine = and) or any one (= or)
    bool require_all = true;

    // Region the tags are restricted to, as placed at the current time; not ok() when there is none
    amrex::RealBox real_box;

    // Height a cell center has when it is tested against the box's z extent
    enum struct BoxZ {
        Grid,        // on the undeformed grid: a static nest, possibly given by indices
        Physical,    // z_phys_cc
        AboveGround  // z_phys_cc less the height of the bottom of its column
    };

    // For a criterion with conditions, the box is only a restriction: it is read from
    //    in_box_lo/in_box_hi into base_box (where it is before any move_* motion), and its z
    //    extent is compared with the physical height or, by default, the height above ground
    //    (erf.<name>.box_z_coordinate). A box given with two values leaves z unrestricted.
    //    Without conditions the box is a static nest and is compared on the undeformed grid.
    BoxZ box_z = BoxZ::Grid;
    amrex::RealBox base_box;

    // Tag only at levels below max_level, and only for min_time <= time <= max_time
    int max_level = 1000;
    amrex::Real min_time = std::numeric_limits<amrex::Real>::lowest();
    amrex::Real max_time = std::numeric_limits<amrex::Real>::max();

    [[nodiscard]] bool is_active (int lev, amrex::Real time) const {
        return (lev >= 0) && (lev < max_level) && (time >= min_time) && (time <= max_time);
    }
};

#endif
