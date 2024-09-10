#ifndef _MYNN_STRUCT_H_
#define _MYNN_STRUCT_H_

enum struct MYNNConfigType {
    NN09 = 0, // Nakanishi & Niino 2009, J Met Soc Japan
    CHEN2021  // Chen & Bryan 2021, JAS
};

struct MYNNLevel25 {

    /*
     * Calculate the stability functions that determine the eddy diffusivities
     * of momentum, heat, QKE, and (optionally) moisture.
     */
    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void calc_stability_funcs(amrex::Real& SM,
                              amrex::Real& SH,
                              amrex::Real& SQ,
                              const amrex::Real GM,
                              const amrex::Real GH,
                              const amrex::Real alphac = 1.0) const
    {
        amrex::Real alphac2 = alphac * alphac;

        // Compute non-dimensional parameters (notation follows NN09)
        amrex::Real Phi1 = 1.0  - alphac2*3.0*A2*B2*(1-C3)*GH;
        amrex::Real Phi2 = 1.0  - alphac2*9.0*A1*A2*(1-C2)*GH;
        amrex::Real Phi3 = Phi1 + alphac2*9.0*A2*A2*(1-C2)*(1-C5)*GH;
        amrex::Real Phi4 = Phi1 - alphac2*12.0*A1*A2*(1-C2)*GH;
        amrex::Real Phi5 = 6.0*alphac*A1*A1*GM;

        // Compute stability functions
        amrex::Real D = Phi2*Phi4 + Phi5*Phi3;
        SM = alphac * A1 * (Phi3 - 3*C1*Phi4) / D;
        SH = alphac * A2 * (Phi2 + 3*C1*Phi5) / D;
        SQ = 3.0 * SM; // revised in NN09
    }

    // Closure coefficients (from Nakanishi & Niino 2009 [NN09])
    amrex::Real A1 = 1.18;
    amrex::Real A2 = 0.665;
    amrex::Real B1 = 24.0;
    amrex::Real B2 = 15.0;
    amrex::Real C1 = 0.137;
    amrex::Real C2 = 0.75;
    amrex::Real C3 = 0.352;
    amrex::Real C4 = 0.0;
    amrex::Real C5 = 0.2;

    // Clipping (defaults from WRF v4.5.0 onward)
    amrex::Real SMmin = 0.0;
    amrex::Real SMmax = 1.0e34;
    amrex::Real SHmin = 0.0;
    amrex::Real SHmax = 4.0;

    MYNNConfigType config;

    // Calculate eddy diffusivity for moisture variables
    bool diffuse_moistvars = false;
};

struct MYNNLevel2 {
    void init_coeffs (amrex::Real A1_lvl25,
                      amrex::Real A2_lvl25,
                      amrex::Real B1,
                      amrex::Real B2,
                      amrex::Real C1,
                      amrex::Real C2,
                      amrex::Real C3,
                      amrex::Real /*C4*/,
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
    amrex::Real calc_Rf (const amrex::Real GM, const amrex::Real GH) const
    {
        AMREX_ASSERT(initialized);
        amrex::Real Ri = -GH/(GM + eps);
        return Ri1*(Ri + Ri2 - std::sqrt(Ri*Ri - Ri3*Ri + Ri2*Ri2));
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real calc_SH (const amrex::Real Rf) const
    {
        AMREX_ASSERT(initialized);
        return 3.0*A2*(gam1+gam2)*(Rfc-Rf)/(1.0-Rf);
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real calc_SM(const amrex::Real Rf) const
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
