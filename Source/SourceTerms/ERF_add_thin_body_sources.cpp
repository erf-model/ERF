#include <AMReX_MultiFab.H>
#include <AMReX_ArrayLim.H>
#include <AMReX_BCRec.H>
#include <AMReX_TableData.H>
#include <AMReX_GpuContainers.H>

#include <ERF_NumericalDiffusion.H>
#include <ERF_PlaneAverage.H>
#include <ERF_TI_slow_headers.H>
#include <ERF_Src_headers.H>
#include <ERF_Utils.H>

using namespace amrex;

/**
 * Function for computing the slow RHS for the evolution equations for the density, potential temperature and momentum.
 *
 * @param[in] xmom_src source terms for x-momentum
 * @param[in] ymom_src source terms for y-momentum
 * @param[in] zmom_src source terms for z-momentum
 * @param[in] xflux_imask_lev thin-body mask on x-faces
 * @param[in] yflux_imask_lev thin-body mask on y-faces
 * @param[in] zflux_imask_lev thin-body mask on z-faces
 * @param[in] thin_xforce_lev x-component of forces on thin immersed bodies
 * @param[in] thin_yforce_lev y-component of forces on thin immersed bodies
 * @param[in] thin_zforce_lev z-component of forces on thin immersed bodies
 */

void add_thin_body_sources ( MultiFab & xmom_src,
                             MultiFab & ymom_src,
                             MultiFab & zmom_src,
                             std::unique_ptr<iMultiFab>& xflux_imask_lev,
                             std::unique_ptr<iMultiFab>& yflux_imask_lev,
                             std::unique_ptr<iMultiFab>& zflux_imask_lev,
                             std::unique_ptr<MultiFab>& thin_xforce_lev,
                             std::unique_ptr<MultiFab>& thin_yforce_lev,
                             std::unique_ptr<MultiFab>& thin_zforce_lev)
{
    BL_PROFILE_REGION("erf_add_thin_body_sources()");

    const bool l_have_thin_xforce = (thin_xforce_lev != nullptr);
    const bool l_have_thin_yforce = (thin_yforce_lev != nullptr);
    const bool l_have_thin_zforce = (thin_zforce_lev != nullptr);

    // *****************************************************************************
    // If a thin immersed body is present, add forcing terms
    // *****************************************************************************
    if (l_have_thin_xforce) {
        MultiFab::Copy(*thin_xforce_lev, xmom_src, 0, 0, 1, 0);
        thin_xforce_lev->mult(-1., 0, 1, 0);
        ApplyInvertedMask(*thin_xforce_lev, *xflux_imask_lev, 0);
        MultiFab::Add(xmom_src, *thin_xforce_lev, 0, 0, 1, 0);
    }

    if (l_have_thin_yforce) {
        MultiFab::Copy(*thin_yforce_lev, ymom_src, 0, 0, 1, 0);
        thin_yforce_lev->mult(-1., 0, 1, 0);
        ApplyInvertedMask(*thin_yforce_lev, *yflux_imask_lev, 0);
        MultiFab::Add(ymom_src, *thin_yforce_lev, 0, 0, 1, 0);
    }

    if (l_have_thin_zforce) {
        MultiFab::Copy(*thin_zforce_lev, zmom_src, 0, 0, 1, 0);
        thin_zforce_lev->mult(-1., 0, 1, 0);
        ApplyInvertedMask(*thin_zforce_lev, *zflux_imask_lev, 0);
        MultiFab::Add(zmom_src, *thin_zforce_lev, 0, 0, 1, 0);
    }

#if 0
#ifndef AMREX_USE_GPU
    if (l_have_thin_xforce) {
        // TODO: Implement particles to better track and output these data
        if (nrk==2) {
            for ( MFIter mfi(S_data[IntVars::cons],TileNoZ()); mfi.isValid(); ++mfi)
            {
                const Box& tbx = mfi.nodaltilebox(0);
                const Array4<const Real> & fx = thin_xforce_lev->const_array(mfi);
                const Array4<const int> & mask = xflux_imask_lev->const_array(mfi);
                ParallelFor(tbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    if (mask(i,j,k)==0) {
                        amrex::AllPrint() << "thin body fx"<<IntVect(i,j,k)<<" = " << fx(i,j,k) << std::endl;
                    }
                });
            }
        }
    }
#endif
#endif

#if 0
#ifndef AMREX_USE_GPU
    if (l_have_thin_yforce) {
        // TODO: Implement particles to better track and output these data
        if (nrk==2) {
            for ( MFIter mfi(S_data[IntVars::cons],TileNoZ()); mfi.isValid(); ++mfi)
            {
                const Box& tby = mfi.nodaltilebox(1);
                const Array4<const Real> & fy = thin_yforce_lev->const_array(mfi);
                const Array4<const int> & mask = yflux_imask_lev->const_array(mfi);
                ParallelFor(tby, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    if (mask(i,j,k)==0) {
                        amrex::AllPrint() << "thin body fy"<<IntVect(i,j,k)<<" = " << fy(i,j,k) << std::endl;
                    }
                });
            }
        }
    }
#endif
#endif

#if 0
#ifndef AMREX_USE_GPU
    if (l_have_thin_zforce) {
        // TODO: Implement particles to better track and output these data
        if (nrk==2) {
            for ( MFIter mfi(S_data[IntVars::cons],TileNoZ()); mfi.isValid(); ++mfi)
            {
                const Box& tbz = mfi.nodaltilebox(2);
                const Array4<const Real> & fz = thin_zforce_lev->const_array(mfi);
                const Array4<const int> & mask = zflux_imask_lev->const_array(mfi);
                ParallelFor(tbz, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    if (mask(i,j,k)==0) {
                        amrex::AllPrint() << "thin body fz"<<IntVect(i,j,k)<<" = " << fz(i,j,k) << std::endl;
                    }
                });
            }
        }
    }
#endif
#endif
}
