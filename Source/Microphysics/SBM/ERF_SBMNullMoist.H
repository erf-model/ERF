#ifndef ERF_SBM_NULL_MOIST_H_
#define ERF_SBM_NULL_MOIST_H_

#include "ERF_NullMoist.H"

/** No-op P1 microphysics adapter; transport is owned by the ERF auxiliary hook. */
class SBMNullMoist final : public NullMoist {
public:
    int Qstate_Moist_Size() override { return 3; } // qv, qc, qr in cons
    int Qstate_Moist_NumConc_Size() override { return 0; }
    int Qstate_NonMoist_Size() override { return 0; }
};

#endif
