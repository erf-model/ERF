#include <ERF_ParFunctions.H>

/**
 * Initialize a Rayleigh damping layer with the same structure as in WRF, based
 * on Durran and Klemp 1983
*/
void
erf_init_rayleigh (amrex::Vector<amrex::Vector<amrex::Real> >& rayleigh_ptrs,
                   amrex::Geometry const& geom,
                   std::unique_ptr<amrex::MultiFab>& /*z_phys_nd*/,
                   amrex::Real /*zdamp*/) override
{
    amrex::ParmParse pp("prob");
    amrex::Real U_0 =   0.0; pp.query("rayleigh_U_0", U_0);
    amrex::Real V_0 =   0.0; pp.query("rayleigh_V_0", V_0);
    amrex::Real W_0 =   0.0; pp.query("rayleigh_W_0", W_0);
    amrex::Real T_0 = 300.0; pp.query("rayleigh_T_0", T_0);

    const int khi = geom.Domain().bigEnd()[2];

    for (int k = 0; k <= khi; k++)
    {
        rayleigh_ptrs[Rayleigh::ubar][k]     = U_0;
        rayleigh_ptrs[Rayleigh::vbar][k]     = V_0;
        rayleigh_ptrs[Rayleigh::wbar][k]     = W_0;
        rayleigh_ptrs[Rayleigh::thetabar][k] = T_0;
    }
}
