#ifndef ERF_SBM_NULL_MOIST_H_
#define ERF_SBM_NULL_MOIST_H_

#include "ERF_NullMoist.H"

/** No-op adapter: ERF keeps qv while SBM owns the projected qc/qr state. */
class SBMNullMoist final : public NullMoist {
public:
    int Qstate_Moist_Size() override { return 3; }
    int Qstate_Moist_NumConc_Size() override { return 0; }
};

#endif
