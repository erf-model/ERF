/**
 * Initialize a Rayleigh damping layer with the same structure as in WRF, based
 * on Durran and Klemp 1983
*/
void
erf_init_rayleigh (amrex::Vector<amrex::Vector<amrex::Real> >& rayleigh_ptrs,
                   amrex::Geometry      const& geom,
                   std::unique_ptr<amrex::MultiFab>& z_phys_cc) override
{
    const int khi = geom.Domain().bigEnd()[2];
    amrex::Vector<amrex::Real> zcc(khi+1);

    if (z_phys_cc) {
        // use_terrain=1
        // calculate the damping strength based on the max height at each k
        amrex::MultiArray4<const amrex::Real> const& ma_z_phys = z_phys_cc->const_arrays();
        for (int k = 0; k <= khi; k++) {
            zcc[k] = amrex::ParReduce(amrex::TypeList<amrex::ReduceOpMax>{},
                                      amrex::TypeList<amrex::Real>{},
                                      *z_phys_cc,
                                      amrex::IntVect(0),
            [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int) noexcept
                -> amrex::GpuTuple<amrex::Real>
            {
                return { ma_z_phys[box_no](i,j,k) };
            });
        }
        amrex::ParallelDescriptor::ReduceRealMax(zcc.data(), zcc.size());

    } else {
        const amrex::Real* prob_lo = geom.ProbLo();
        const auto dx              = geom.CellSize();
        for (int k = 0; k <= khi; k++)
        {
            zcc[k] = prob_lo[2] + (k+0.5) * dx[2];
        }
    }

    const auto ztop = zcc[khi];
    amrex::Print() << "Rayleigh damping (coef="<<parms.dampcoef<<" s^-1)"
        << " between " << ztop-parms.zdamp << " and " << ztop << " m"
        << std::endl;

    for (int k = 0; k <= khi; k++)
    {
        const amrex::Real zfrac = 1 - (ztop - zcc[k]) / parms.zdamp;
        if (zfrac >= 0)
        {
            const amrex::Real sinefac = std::sin(PIoTwo*zfrac);
            rayleigh_ptrs[Rayleigh::tau][k]      = parms.dampcoef * sinefac * sinefac;
            rayleigh_ptrs[Rayleigh::ubar][k]     = parms.U_0;
            rayleigh_ptrs[Rayleigh::vbar][k]     = parms.V_0;
            rayleigh_ptrs[Rayleigh::wbar][k]     = parms.W_0;
            rayleigh_ptrs[Rayleigh::thetabar][k] = parms.T_0;
        }
        else
        {
            rayleigh_ptrs[Rayleigh::tau][k]      = amrex::Real(0.0);
            rayleigh_ptrs[Rayleigh::ubar][k]     = amrex::Real(0.0);
            rayleigh_ptrs[Rayleigh::vbar][k]     = amrex::Real(0.0);
            rayleigh_ptrs[Rayleigh::wbar][k]     = amrex::Real(0.0);
            rayleigh_ptrs[Rayleigh::thetabar][k] = amrex::Real(0.0);
        }
    }
}
