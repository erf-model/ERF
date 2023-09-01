/* Initialize a Rayleigh damping layer with the same structure as in WRF, based
 * on Durran and Klemp 1983
 */
void
erf_init_rayleigh(Vector<Real>& tau,
                  Vector<Real>& ubar,
                  Vector<Real>& vbar,
                  Vector<Real>& wbar,
                  Vector<Real>& thetabar,
                  Geometry      const& geom)
{
    const auto ztop = geom.ProbHi()[2];
    amrex::Print() << "Rayleigh damping (coef="<<parms.dampcoef<<") between "
      << ztop-parms.zdamp << " and " << ztop << std::endl;
    
    const Real* prob_lo = geom.ProbLo();
    const auto dx       = geom.CellSize();
    const int khi       = geom.Domain().bigEnd()[2];
    
    for (int k = 0; k <= khi; k++)
    {
        const Real z = prob_lo[2] + (k + 0.5) * dx[2];
        const Real zfrac = 1 - (ztop - z) / parms.zdamp;
        if (zfrac >= 0)
        {
            const Real sinefac = std::sin(PIoTwo*zfrac);
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
