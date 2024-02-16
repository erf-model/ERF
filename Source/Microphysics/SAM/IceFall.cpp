#include <AMReX_ParReduce.H>
#include "SAM.H"
#include "TileNoZ.H"

using namespace amrex;

/**
 * Computes contributions to Microphysics and thermodynamic variables from falling cloud ice in each column.
 */
void SAM::IceFall () {

    Real dz  = m_geom.CellSize(2);
    Real dtn = dt;
    int  nz  = nlev;

    int kmax, kmin;
    auto qcl   = mic_fab_vars[MicVar::qcl];
    auto qci   = mic_fab_vars[MicVar::qci];
    auto qt    = mic_fab_vars[MicVar::qt];
    auto tabs  = mic_fab_vars[MicVar::tabs];
    auto theta = mic_fab_vars[MicVar::theta];

    MultiFab fz;
    fz.define(qcl->boxArray(),qcl->DistributionMap(), 1, qcl->nGrowVect());
    fz.setVal(0.);

    auto qifall_t = qifall.table();
    auto tlatqi_t = tlatqi.table();
    auto rho1d_t  = rho1d.table();

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

//std::cout << "ice_fall: " << kmin << "; " << kmax << std::endl;

    // for (int k=0; k<nzm; k++) {
    //  for (int icrm=0; icrm<ncrms; icrm++) {
    ParallelFor(nz, [=] AMREX_GPU_DEVICE (int k)
    {
        qifall_t(k) = 0.0;
        tlatqi_t(k) = 0.0;
    });

    //if(index_cloud_ice == -1) { return;}

    // for (int k=0; k<nzm; k++) {
    //   for (int j=0; j<ny; j++) {
    //     for (int i=0; i<nx; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    for ( amrex::MFIter mfi(*tabs, TileNoZ()); mfi.isValid(); ++mfi) {
        //auto qcl_array   = qcl->array(mfi);
        auto qci_array   = qci->array(mfi);
        auto qt_array    = qt->array(mfi);
        //auto tabs_array  = tabs->array(mfi);
        auto theta_array = theta->array(mfi);
        auto fz_array    = fz.array(mfi);

        const auto& box3d = mfi.tilebox();

        ParallelFor(box3d, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            if (k >= std::max(0,kmin-1) && k <= kmax ) {
                // Set up indices for x-y planes above and below current plane.
                int kc = std::min(k+1, nz-1);
                int kb = std::max(k-1, 0);

                // CFL number based on grid spacing interpolated to interface i,j,k-1/2
                Real coef = dtn/dz; //dtn/(0.5*(adz(kb)+adz(k))*dz);

                // Compute cloud ice density in this cell and the ones above/below.
                // Since cloud ice is falling, the above cell is u(icrm,upwind),
                // this cell is c (center) and the one below is d (downwind).
                Real qiu = rho1d_t(kc)*qci_array(i,j,kc);
                Real qic = rho1d_t(k )*qci_array(i,j,k );
                Real qid = rho1d_t(kb)*qci_array(i,j,kb);

                // Ice sedimentation velocity depends on ice content. The fiting is
                // based on the data by Heymsfield (JAS,2003). -Marat
                Real vt_ice = min( 0.4 , 8.66 * pow( (max(0.,qic)+1.e-10) , 0.24) );   // Heymsfield (JAS, 2003, p.2607)

                // Use MC flux limiter in computation of flux correction.
                // (MC = monotonized centered difference).
                //         if (qic.eq.qid) then
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

        // for (int j=0; j<ny; j++) {
        //  for (int i=0; i<nx; i++) {
        //    for (int icrm=0; icrm<ncrms; icrm++) {
        /* by Xingqiu Yuan, need setup boundary condition for fz_t here???
           ParallelFor(box2d, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
           fz_t(nz-1,j,i) = 0.0;
           });
        */

        // for (int k=0; k<nzm; k++) {
        //   for (int j=0; j<ny; j++) {
        //     for (int i=0; i<nx; i++) {
        //       for (int icrm=0; icrm<ncrms; icrm++) {
        Real fac_cond = m_fac_cond;
        Real fac_fus  = m_fac_fus;
        ParallelFor(box3d, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            if ( k >= std::max(0,kmin-1) && k <= kmax ) {
                Real coef = dtn/dz;
                // The cloud ice increment is the difference of the fluxes.
                Real dqi  = coef*(fz_array(i,j,k)-fz_array(i,j,k+1));
                // Add this increment to both non-precipitating and total water.
                amrex::Gpu::Atomic::Add(&qt_array(i,j,k), dqi);
                // Include this effect in the total moisture budget.
                amrex::Gpu::Atomic::Add(&qifall_t(k), dqi);

                // The latent heat flux induced by the falling cloud ice enters
                // the liquid-ice static energy budget in the same way as the
                // precipitation.  Note: use latent heat of sublimation.
                Real lat_heat = (fac_cond+fac_fus)*dqi;

                // Add divergence of latent heat flux contribution to liquid-ice static potential temperature.
                amrex::Gpu::Atomic::Add(&theta_array(i,j,k), -lat_heat);
                // Add divergence to liquid-ice static energy budget.
                amrex::Gpu::Atomic::Add(&tlatqi_t(k), -lat_heat);
            }
        });

#if 0
        // for (int j=0; j<ny; j++) {
        //    for (int i=0; i<nx; i++) {
        //      for (int icrm=0; icrm<ncrms; icrm++) {
        //parallel_for( SimpleBounds<3>(ny,nx,ncrms) , YAKL_LAMBDA (int j, int i, int icrm) {
        ParallelFor(box2d, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            Real coef = dtn/dz;
            Real dqi = -coef*fz(i,j,0);
            precsfc (i,j) = precsfc (i,j)+dqi;
            precssfc(i,j) = precssfc(i,j)+dqi;
        });
#endif
    }
}

