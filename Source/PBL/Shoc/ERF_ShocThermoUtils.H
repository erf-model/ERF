#ifndef ERF_SHOC_THERMO_UTILS_H_
#define ERF_SHOC_THERMO_UTILS_H_

#include "ERF_Constants.H"
#include "ERF_ShocConstants.H"
#include "ERF_ShocTypes.H"

namespace shoc {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
latent_sublimation () noexcept
{
    return L_v + constants::latent_ice();
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
condensate_latent_heat (amrex::Real qc,
                        amrex::Real qi) noexcept
{
    return L_v * amrex::max(amrex::Real(0.0), qc)
         + latent_sublimation() * amrex::max(amrex::Real(0.0), qi);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
theta_from_thetal (amrex::Real thetal,
                   amrex::Real qc,
                   amrex::Real qi,
                   amrex::Real exner) noexcept
{
    return thetal
         + condensate_latent_heat(qc, qi)
           / (Cp_d * amrex::max(exner, amrex::Real(1.0e-12)));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
temperature_from_thetal (amrex::Real thetal,
                         amrex::Real qc,
                         amrex::Real qi,
                         amrex::Real exner) noexcept
{
    return amrex::max(constants::min_temp(),
                      thetal * amrex::max(exner, amrex::Real(1.0e-12))
                      + condensate_latent_heat(qc, qi) / Cp_d);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
thetal_from_theta (amrex::Real theta,
                   amrex::Real qc,
                   amrex::Real qi,
                   amrex::Real exner) noexcept
{
    return theta
         - condensate_latent_heat(qc, qi)
           / (Cp_d * amrex::max(exner, amrex::Real(1.0e-12)));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
moist_energy (amrex::Real tabs,
              amrex::Real z,
              amrex::Real qv,
              amrex::Real qc,
              amrex::Real qi,
              amrex::Real u,
              amrex::Real v,
              amrex::Real tke) noexcept
{
    // Ice-reference moist energy has no explicit qi term; qi is retained in
    // the signature so callers pass a complete vapor/liquid/ice state.
    amrex::ignore_unused(qi);
    return Cp_d * tabs
         + CONST_GRAV * z
         + amrex::Real(0.5) * (u * u + v * v)
         + tke
         + latent_sublimation() * qv
         + constants::latent_ice() * amrex::max(amrex::Real(0.0), qc);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
reconstruct_pdf_state (amrex::Real thetal,
                       amrex::Real qw,
                       amrex::Real exner,
                       amrex::Real qi_seed,
                       amrex::Real pdf_ql,
                       amrex::Real& tabs,
                       amrex::Real& qv,
                       amrex::Real& qc,
                       amrex::Real& qi) noexcept
{
    // Interim native-SHOC contract: preserve any existing ice seed for
    // phase-aware thermodynamics, but do not create or repartition cloud ice.
    const amrex::Real qi_total = shoc_clamp(qi_seed, amrex::Real(0.0), amrex::max(amrex::Real(0.0), qw));
    const amrex::Real ql_total = shoc_clamp(pdf_ql, amrex::Real(0.0), amrex::max(amrex::Real(0.0), qw - qi_total));

    qc = ql_total;
    qi = qi_total;
    tabs = temperature_from_thetal(thetal, qc, qi, exner);
    qv = amrex::max(amrex::Real(0.0), qw - qc - qi);
}

} // namespace shoc

#endif
