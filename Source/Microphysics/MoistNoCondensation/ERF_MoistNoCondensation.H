#ifndef ERF_MOISTNOCOND_H
#define ERF_MOISTNOCOND_H

/*
 * MoistNoCondensation is a shell moisture model that does not modify anything.
 *
 * The scheme allows qv to be allocated and updated by the dycore without any
 * source terms due to macro- or microphysics.
 */

#include "ERF_NullMoist.H"

class MoistNoCondensation : public NullMoist {

public:
    // constructor
    MoistNoCondensation () {}

    // destructor
    virtual ~MoistNoCondensation () = default;

    int
    Qstate_Moist_Size () override { return MoistNoCondensation::n_qstate_moist_size; }

    int
    Qstate_Moist_NumConc_Size () override { return MoistNoCondensation::n_qstate_moist_numconc_size; }

private:

    // Number of qstate variables
    int n_qstate_moist_size = 2;

    // Number of qstate variables that are number concentrations
    int n_qstate_moist_numconc_size = 0;
};
#endif
