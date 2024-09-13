#include <Utils/ERF_ParFunctions.H>

/**
 * Initialize a Rayleigh damping layer with the same structure as in WRF, based
 * on Durran and Klemp 1983
*/
void
erf_init_rayleigh (amrex::Vector<amrex::Vector<amrex::Real> >& rayleigh_ptrs,
                   amrex::Geometry const& geom,
                   std::unique_ptr<amrex::MultiFab>& z_phys_nd,
                   amrex::Real zdamp) override
{
    const int khi = geom.Domain().bigEnd()[2];
    for (int k = 0; k <= khi; k++)
    {
        rayleigh_ptrs[Rayleigh::ubar][k]     = parms.U_0;
        rayleigh_ptrs[Rayleigh::vbar][k]     = parms.V_0;
        rayleigh_ptrs[Rayleigh::wbar][k]     = parms.W_0;
        rayleigh_ptrs[Rayleigh::thetabar][k] = parms.T_0;
    }
}
