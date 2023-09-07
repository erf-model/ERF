/**
 * Initialize a Rayleigh damping layer with the same structure as in WRF, based
 * on Durran and Klemp 1983
*/
void
erf_init_rayleigh(amrex::Vector<amrex::Real>& tau,
                  amrex::Vector<amrex::Real>& ubar,
                  amrex::Vector<amrex::Real>& vbar,
                  amrex::Vector<amrex::Real>& wbar,
                  amrex::Vector<amrex::Real>& thetabar,
                  amrex::Geometry      const& geom) override
{
    const auto ztop = geom.ProbHi()[2];
    amrex::Print() << "Rayleigh damping (coef="<<parms.dampcoef<<") between "
      << ztop-parms.zdamp << " and " << ztop << std::endl;

    const amrex::Real* prob_lo = geom.ProbLo();
    const auto dx              = geom.CellSize();
    const int khi              = geom.Domain().bigEnd()[2];

    for (int k = 0; k <= khi; k++)
    {
        const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];
        const amrex::Real zfrac = 1 - (ztop - z) / parms.zdamp;
        if (zfrac >= 0)
        {
            const amrex::Real sinefac = std::sin(PIoTwo*zfrac);
            tau[k]      = parms.dampcoef * sinefac * sinefac;
            ubar[k]     = parms.U_0;
            vbar[k]     = parms.V_0;
            wbar[k]     = parms.W_0;
            thetabar[k] = parms.T_0;
        }
        else
        {
            tau[k]      = 0.0;
            ubar[k]     = 0.0;
            vbar[k]     = 0.0;
            wbar[k]     = 0.0;
            thetabar[k] = 0.0;
        }
    }
}
