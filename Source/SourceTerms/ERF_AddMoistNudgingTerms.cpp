#if defined(ERF_USE_NETCDF)

#include <ERF_SrcHeaders.H>

using namespace amrex;

/**
 * Function for computing the slow RHS for the evolution equations for the density, potential temperature and momentum.
 *
 * @param[in] S_data  current solution
 * @param[in] source  source terms for moist conserved variables
 * @param[in] dt      slow time step
 * @param[in] old_stage_time_total
 * @param[in] start_bdy_time
 * @param[in] final_bdy_time
 * @param[in] bdy_time_interval
 * @param[in] bdy_factor
 * @param[in] width
 * @param[in] geom   Container for geometric information
 * @param[in] bdy_data_xlo boundary data on interior of low x-face
 * @param[in] bdy_data_xhi boundary data on interior of high x-face
 * @param[in] bdy_data_ylo boundary data on interior of low y-face
 * @param[in] bdy_data_yhi boundary data on interior of high y-face
 * @param[in] m_r2d        boundary data read in using ReadBndryPlane
 */

void add_moist_nudging_terms (const MultiFab& S_data,
                                    MultiFab & source,
                              const Real& dt,
                              const Real& old_stage_time_total,
                              const Real& start_bdy_time,
                              const Real& final_bdy_time,
                              const Real& bdy_time_interval,
                              const Real& bdy_factor,
                              int  width,
                              const Geometry& geom,
                              Vector<Vector<FArrayBox>>& bdy_data_xlo,
                              Vector<Vector<FArrayBox>>& bdy_data_xhi,
                              Vector<Vector<FArrayBox>>& bdy_data_ylo,
                              Vector<Vector<FArrayBox>>& bdy_data_yhi,
                              std::unique_ptr<ReadBndryPlanes>& m_r2d)
{
    BL_PROFILE_REGION("erf_add_moist_nudging_terms()");

    const Box domain = geom.Domain();

    // Temporary MF so we can nudge qv + qc to the bdy data
    int nc = S_data.nComp();
    IntVect ng  = S_data.nGrowVect();
    BoxArray ba = S_data.boxArray();
    DistributionMapping dm = S_data.DistributionMap();
    MultiFab S_tmp(ba, dm, nc, ng);
    MultiFab::Copy(S_tmp, S_data, Rho_comp  , Rho_comp  , 1, ng);
    MultiFab::Copy(S_tmp, S_data, RhoQ1_comp, RhoQ1_comp, 1, ng);
    MultiFab::Add (S_tmp, S_data, RhoQ2_comp, RhoQ1_comp, 1, ng);
    if (nc > RhoQ6_comp) {
        MultiFab::Add (S_tmp, S_data, RhoQ3_comp, RhoQ1_comp, 1, ng);
    }

    for (MFIter mfi(S_data,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box tbx  = mfi.tilebox();

        const Array4<const Real>& new_cons_const = S_tmp.const_array(mfi);
        const Array4<      Real>& src_arr        = source.array(mfi);
        //
        // Note that old_stage_time_total = start_time+old_stage_time is total time
        //           start_bdy_time and final_bdy_time are total time
        //
        moist_set_rhs(geom, tbx, new_cons_const, src_arr,
                      old_stage_time_total, dt,
                      start_bdy_time, final_bdy_time, bdy_time_interval,
                      bdy_factor, width, domain,
                      bdy_data_xlo, bdy_data_xhi,
                      bdy_data_ylo, bdy_data_yhi,
                      m_r2d);
    }
}
#endif
