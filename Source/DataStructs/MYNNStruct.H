#ifndef _MYNN_STRUCT_H_
#define _MYNN_STRUCT_H_

struct MYNNLevel2 {
    void init_coeffs (amrex::Real A1_lvl25,
                      amrex::Real A2_lvl25,
                      amrex::Real B1,
                      amrex::Real B2,
                      amrex::Real C1,
                      amrex::Real C2,
                      amrex::Real C3,
                      amrex::Real C4,
                      amrex::Real C5)
    {
        A1 = A1_lvl25;
        A2 = A2_lvl25;
        gam2 = (2.0*A1*(3.0-2.0*C2) + B2*(1.0-C3)) / B1;
        F1  = B1*(gam1-C1) + 2.0*A1*(3.0-2.0*C2) + 3.0*A2*(1.0-C2)*(1.0-C5);
        F2  = B1*(gam1+gam2) - 3.0*A1*(1.0-C2);
        Rf1 = B1*(gam1-C1) / F1;
        Rf2 = B1*gam1 / F2;
        Rfc = gam1 / (gam1 + gam2);
        Ri1 = 0.5*A2*F2/(A1*F1);
        Ri2 = 0.5*Rf1/Ri1;
        Ri3 = (2.0*Rf2-Rf1)/Ri1;

        initialized = true;
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real calc_Rf (amrex::Real GM, amrex::Real GH) const
    {
        AMREX_ASSERT(initialized);
        amrex::Real Ri = -GH/(GM + eps);
        return Ri1*(Ri + Ri2 - std::sqrt(Ri*Ri - Ri3*Ri + Ri2*Ri2));
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real calc_SH (amrex::Real Rf) const
    {
        AMREX_ASSERT(initialized);
        return 3.0*A2*(gam1+gam2)*(Rfc-Rf)/(1.0-Rf);
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real calc_SM(amrex::Real Rf) const
    {
        AMREX_ASSERT(initialized);
        return A1*F1/(A2*F2) * (Rf1-Rf) / (Rf2-Rf) * calc_SH(Rf);
    }

    bool initialized{false};

    amrex::Real A1, A2; // from Level 2.5
    amrex::Real F1, F2;
    const amrex::Real gam1 = 0.235;
    amrex::Real gam2;
    amrex::Real Rf1, Rf2, Rfc;
    amrex::Real Ri1, Ri2, Ri3;

    const amrex::Real eps = std::numeric_limits<amrex::Real>::epsilon();
};
#endif
