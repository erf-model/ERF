#ifndef ERF_MYNN_STRUCT_H_
#define ERF_MYNN_STRUCT_H_

#include <AMReX_Enum.H>

#include "ERF_Constants.H"

// NN09     : Nakanishi & Niino 2009, J Met Soc Japan -- Lt_alpha = 0.23
// Chen2021 : Chen & Bryan 2021, JAS                  -- Lt_alpha = 0.10, and
//            the master length scale is the inverse-square (rather than
//            inverse) sum of l_S, l_T and l_B.
// NN09 must stay first so that the zero value remains the default.
AMREX_ENUM(MYNNConfigType, NN09, Chen2021);

/**
 * @brief Coefficients and stability functions for the Level 2.5 MYNN closure.
 */
struct MYNNLevel25 {

    /*
     * Calculate the stability functions that determine the eddy diffusivities
     * of momentum, heat, KE, and (optionally) moisture.
     */
    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void calc_stability_funcs (amrex::Real& SM,
                               amrex::Real& SH,
                               amrex::Real& SQ,
                               const amrex::Real GM,
                               const amrex::Real GH,
                               const amrex::Real alphac = one) const
    {
        amrex::Real alphac2 = alphac * alphac;

        // Compute non-dimensional parameters (notation follows NN09, Eqns. 33-37)
        amrex::Real Phi1 = one  -  three * alphac2 * A2 * B2 * (1-C3) * GH;
        amrex::Real Phi2 = one  -  amrex::Real(9.0) * alphac2 * A1 * A2 * (1-C2) * GH;
        amrex::Real Phi3 = Phi1 +  amrex::Real(9.0) * alphac2 * A2 * A2 * (1-C2) * (1-C5) * GH;
        amrex::Real Phi4 = Phi1 - amrex::Real(12.0) * alphac2 * A1 * A2 * (1-C2) * GH;
        amrex::Real Phi5 =         amrex::Real(6.0) * alphac2 * A1 * A1 * GM;

        // Compute level amrex::Real(2.5) stability functions
        amrex::Real D = Phi2*Phi4 + Phi5*Phi3;      // NN09 Eqn. 31
        SM = alphac * A1 * (Phi3 - 3*C1*Phi4) / D;  // NN09 Eqn. 27
        SH = alphac * A2 * (Phi2 + 3*C1*Phi5) / D;  // NN09 Eqn. 28
        SQ = SQfac  * SM;                           // NN09 Eqn. 67
    }

    // Closure coefficients (from Nakanishi & Niino 2009 [NN09])
    amrex::Real A1 = amrex::Real(1.18);
    amrex::Real A2 = amrex::Real(0.665);
    amrex::Real B1 = amrex::Real(24.0);
    amrex::Real B2 = amrex::Real(15.0);
    amrex::Real C1 = amrex::Real(0.137);
    amrex::Real C2 = amrex::Real(0.75);
    amrex::Real C3 = amrex::Real(0.352);
    amrex::Real C4 = zero;
    amrex::Real C5 = amrex::Real(0.2);

    // Third-order turbulent fluxes
    amrex::Real SQfac = three; // revised from two in NN09

    // Clipping (defaults from WRF v4.5.0 onward)
    amrex::Real SMmin = zero;
    amrex::Real SMmax = bogus_large_value;
    amrex::Real SHmin = zero;
    amrex::Real SHmax = amrex::Real(4);
    amrex::Real SQmin = zero;
    amrex::Real SQmax = SQfac * SHmax;

    MYNNConfigType config = MYNNConfigType::NN09;

    // Coefficient on the ABL-depth length scale, l_T = Lt_alpha*<zq>/<q>
    // (NN09, Eqn. 54).  Set from `config` in TurbChoice::init_params -- 0.23
    // for NN09, 0.10 for Chen2021 -- and then overridable on its own with
    // erf.pbl_mynn_Lt_alpha.  Mellor-Yamada 1982, which COAMPS uses, takes 0.1.
    amrex::Real Lt_alpha = amrex::Real(0.23);

    // Optional taper of the ABL-depth length scale toward the top of the
    // boundary layer: l_T -> l_T * (1 - z/h)^Lt_taper_exp, with h the
    // diagnosed PBL height.  Zero (the default) disables it and recovers the
    // untapered NN09 form exactly.
    //
    // l_T (NN09 Eqn. 54) is a single number for the whole column, so once
    // kappa*z exceeds it the master length scale saturates and stays flat all
    // the way to the inversion.  l_S only grows with height and l_B only bites
    // once the layer is stably stratified, so in a mixed layer that is neutral
    // or slightly unstable nothing brings the length scale back down.
    amrex::Real Lt_taper_exp = amrex::Real(0.0);
};

/**
 * @brief Coefficients and stability functions for the Level 2 MYNN closure.
 */
struct MYNNLevel2 {
    /**
     * @brief Initialize Level 2 closure coefficients from Level 2.5 coefficients.
     * @param[in] A1_lvl25 Level 2.5 coefficient A1.
     * @param[in] A2_lvl25 Level 2.5 coefficient A2.
     * @param[in] B1 Closure coefficient B1.
     * @param[in] B2 Closure coefficient B2.
     * @param[in] C1 Closure coefficient C1.
     * @param[in] C2 Closure coefficient C2.
     * @param[in] C3 Closure coefficient C3.
     * @param[in] C5 Closure coefficient C5.
     */
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
        gam2 = (two*A1*(three-two*C2) + B2*(one-C3)) / B1;
        F1  = B1*(gam1-C1) + two*A1*(three-two*C2) + three*A2*(one-C2)*(one-C5);
        F2  = B1*(gam1+gam2) - three*A1*(one-C2);
        Rf1 = B1*(gam1-C1) / F1;
        Rf2 = B1*gam1 / F2;
        Rfc = gam1 / (gam1 + gam2);
        Ri1 = myhalf*A2*F2/(A1*F1);
        Ri2 = myhalf*Rf1/Ri1;
        Ri3 = (two*Rf2-Rf1)/Ri1;

        initialized = true;
    }

    /**
     * @brief Calculate the reduced Richardson number.
     * @param[in] GM Momentum gradient.
     * @param[in] GH Heat gradient.
     * @return Reduced Richardson number.
     */
    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real calc_Rf (const amrex::Real GM, const amrex::Real GH) const
    {
        AMREX_ASSERT(initialized);
        amrex::Real lGM = std::copysign(std::max(std::fabs(GM),eps),GM);
        amrex::Real Ri  = -GH/lGM;
        return Ri1*(Ri + Ri2 - std::sqrt(Ri*Ri - Ri3*Ri + Ri2*Ri2));
    }

    /**
     * @brief Calculate the heat stability function.
     * @param[in] Rf Reduced Richardson number.
     * @return Heat stability function.
     */
    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real calc_SH (const amrex::Real Rf) const
    {
        AMREX_ASSERT(initialized);
        return three*A2*(gam1+gam2)*(Rfc-Rf)/(one-Rf);
    }

    /**
     * @brief Calculate the momentum stability function.
     * @param[in] Rf Reduced Richardson number.
     * @return Momentum stability function.
     */
    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real calc_SM (const amrex::Real Rf) const
    {
        AMREX_ASSERT(initialized);
        return A1*F1/(A2*F2) * (Rf1-Rf) / (Rf2-Rf) * calc_SH(Rf);
    }

    bool initialized{false};

    amrex::Real A1, A2; // from Level amrex::Real(2.5)
    amrex::Real F1, F2;
    const amrex::Real gam1 = amrex::Real(0.235);
    amrex::Real gam2;
    amrex::Real Rf1, Rf2, Rfc;
    amrex::Real Ri1, Ri2, Ri3;

    const amrex::Real eps = std::numeric_limits<amrex::Real>::epsilon();
};
#endif
