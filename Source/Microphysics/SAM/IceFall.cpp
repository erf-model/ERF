#include <AMReX_ParReduce.H>
#include "SAM.H"
#include "TileNoZ.H"

using namespace amrex;

/**
 * Sedimentation of cloud ice (A32/33); sources enthalpy
 */
void SAM::IceFall () {

    Real dz  = m_geom.CellSize(2);
    Real dtn = dt;
    int  nz  = nlev;

    int kmax, kmin;
    auto qcl   = mic_fab_vars[MicVar::qcl];
    auto qci   = mic_fab_vars[MicVar::qci];
    auto qn    = mic_fab_vars[MicVar::qn];
    auto qt    = mic_fab_vars[MicVar::qt];
    auto rho   = mic_fab_vars[MicVar::rho];
    auto tabs  = mic_fab_vars[MicVar::tabs];

    MultiFab fz;
    IntVect  ng = qcl->nGrowVect();
    BoxArray ba = qcl->boxArray();
    DistributionMapping dm = qcl->DistributionMap();
    fz.define(convert(ba, IntVect(0,0,1)), dm, 1, ng);
    fz.setVal(0.);

    kmin = nz;
    kmax = -1;

    { // calculate maximum and minium ice fall vertical region
        auto const& qcl_arrays  = qcl->const_arrays();
        auto const& qci_arrays  = qci->const_arrays();
        auto const& tabs_arrays = tabs->const_arrays();

        GpuTuple<int, int> k_max_min = ParReduce(TypeList<ReduceOpMax, ReduceOpMin>{},
                                                 TypeList<int, int>{},
                                                 *qcl, IntVect::TheZeroVector(),
                                                 [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
                                                 -> GpuTuple<int, int>
                                                 {
                                                     bool is_positive       = qcl_arrays[box_no](i,j,k)+qci_arrays[box_no](i,j,k) > 0.0;
                                                     bool smaller_than_zero = tabs_arrays[box_no](i,j,k) < 273.15;
                                                     int mkmin = is_positive && smaller_than_zero ? k : nz-1;
                                                     int mkmax = is_positive && smaller_than_zero ? k : 0;
                                                     return {mkmax, mkmin};
                                                 });
        kmax = amrex::get<0>(k_max_min);
        kmin = amrex::get<1>(k_max_min);
    }

    for ( amrex::MFIter mfi(*tabs, TileNoZ()); mfi.isValid(); ++mfi) {
        auto qci_array   = qci->array(mfi);
        auto qn_array    = qn->array(mfi);
        auto qt_array    = qt->array(mfi);
        auto rho_array   = rho->array(mfi);
        auto fz_array    = fz.array(mfi);

        const auto& gbox3d = mfi.tilebox(IntVect(0),IntVect(0,0,1));
        const auto& box3d  = mfi.tilebox();

        ParallelFor(gbox3d, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            if (k >= std::max(0,kmin-1) && k <= kmax+1 ) {
                // Set up indices for x-y planes above and below current plane.
                int kc = std::min(k+1, nz-1);
                int kb = std::max(k-1, 0);

                // CFL number based on grid spacing interpolated to interface i,j,k-1/2
                Real coef = dtn/dz;

                // Compute cloud ice density in this cell and the ones above/below.
                // Since cloud ice is falling, the above cell is u(icrm,upwind),
                // this cell is c (center) and the one below is d (downwind).
                Real qiu = rho_array(i,j,kc)*qci_array(i,j,kc);
                Real qic = rho_array(i,j,k )*qci_array(i,j,k );
                Real qid = rho_array(i,j,kb)*qci_array(i,j,kb);

                // Ice sedimentation velocity depends on ice content. The fiting is
                // based on the data by Heymsfield (JAS,2003). -Marat
                Real vt_ice = min( 0.4 , 8.66 * pow( (max(0.,qic)+1.e-10) , 0.24) );   // Heymsfield (JAS, 2003, p.2607)

                // Use MC flux limiter in computation of flux correction.
                // (MC = monotonized centered difference).
                Real tmp_phi;
                if ( std::abs(qic-qid) < 1.0e-25 ) {  // when qic, and qid is very small, qic_qid can still be zero
                    // even if qic is not equal to qid. so add a fix here +++mhwang
                    tmp_phi = 0.;
                } else {
                    Real tmp_theta = (qiu-qic)/(qic-qid+1.0e-20);
                    tmp_phi = max(0., min(0.5*(1.+tmp_theta), min(2., 2.*tmp_theta)));
                }

                // Compute limited flux.
                // Since falling cloud ice is a 1D advection problem, this
                // flux-limited advection scheme is monotonic.
                fz_array(i,j,k) = -vt_ice*(qic - 0.5*(1.-coef*vt_ice)*tmp_phi*(qic-qid));
            }
        });

        ParallelFor(box3d, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            if ( k >= std::max(0,kmin) && k <= kmax ) {
                //==================================================
                // Cloud ice sedimentation (A32)
                //==================================================
                Real coef = dtn/dz;
                Real dqi  = std::max(-qci_array(i,j,k),coef*(fz_array(i,j,k)-fz_array(i,j,k+1)));

                // Add this increment to both non-precipitating and total water.
                qci_array(i,j,k) += dqi;
                 qn_array(i,j,k) += dqi;
                 qt_array(i,j,k) += dqi;

                // NOTE: Sedimentation does not affect the potential temperature,
                //       but it does affect the liquid/ice static energy.
                //       No source to Theta occurs here.
            }
        });
    }
}

