/**
 * \file ERF_FillZeroAreaFaceFluxes.cpp
 */
#include "ERF.H"
#include "ERF_Utils.H"

using namespace amrex;

/**
 * Compute phi gradients where the area of the face is zero
 *
 * At such a face one of the two adjacent cells is covered, so the gradient is extrapolated
 * from the two nearest face gradients on the uncovered side, which reaches three cells away.
 * Give phi three filled ghost cells to get that stencil everywhere; with fewer, faces near
 * the edge of a box fall back to a first-order one-sided gradient, and the result then
 * depends on the grid decomposition.
 *
 * @tparam T EB factory or auxiliary EB data type for face-centered grids
 * @param phi Cell-centered solution used to compute gradients, with its ghost cells filled
 * @param fluxes Face-centered gradient fluxes to fill
 * @param geom Geometry used for inverse cell spacing
 * @param ebfact Cell-centered embedded-boundary factory
 * @param ebfact_u Embedded-boundary data on x-faces
 * @param ebfact_v Embedded-boundary data on y-faces
 * @param ebfact_w Embedded-boundary data on z-faces
 */
template <typename T>
void
FillZeroAreaFaceFluxes (MultiFab& phi, Array<MultiFab,AMREX_SPACEDIM>& fluxes,
                        const Geometry& geom, EBFArrayBoxFactory const& ebfact,
                        T const& ebfact_u, T const& ebfact_v, T const& ebfact_w)
{
    BL_PROFILE("ERF::FillZeroAreaFaceFluxes()");

    const GpuArray<Real, AMREX_SPACEDIM> dxInv = geom.InvCellSizeArray();

    for (MFIter mfi(phi,TileNoZ()); mfi.isValid(); ++mfi)
    {
        const Box& tbx = mfi.tilebox();
        const Box& xbx = mfi.nodaltilebox(0);
        const Box& ybx = mfi.nodaltilebox(1);
        const Box& zbx = mfi.nodaltilebox(2);

        // EBCellFlagFab const& cflag_fab = (ebfact.get_const_factory())->getMultiEBCellFlagFab()[mfi];
        EBCellFlagFab const& cflag_fab = (ebfact.getMultiEBCellFlagFab())[mfi];
        Array4<const EBCellFlag> cflag = cflag_fab.const_array();

        if (cflag_fab.getType(tbx) == FabType::singlevalued)
        {
            Array4<const Real> apx = ebfact.getAreaFrac()[0]->const_array(mfi);
            Array4<const Real> apy = ebfact.getAreaFrac()[1]->const_array(mfi);
            Array4<const Real> apz = ebfact.getAreaFrac()[2]->const_array(mfi);

            Array4<const EBCellFlag> u_cflag = ebfact_u.getMultiEBCellFlagFab()[mfi].const_array();
            Array4<const EBCellFlag> v_cflag = ebfact_v.getMultiEBCellFlagFab()[mfi].const_array();
            Array4<const EBCellFlag> w_cflag = ebfact_w.getMultiEBCellFlagFab()[mfi].const_array();

            Array4<Real const> const& p_arr = phi.const_array(mfi);
            Array4<Real> const& fx = fluxes[0].array(mfi);
            Array4<Real> const& fy = fluxes[1].array(mfi);
            Array4<Real> const& fz = fluxes[2].array(mfi);

            // The region of phi that we may read: this FAB's own cells, less any that lie
            // outside the domain in a non-periodic direction, since the FillBoundary done
            // before this call leaves those unfilled. The three-cell stencils below reach
            // three cells past the tile at a face on the box edge, so they must be checked
            // against this rather than assumed to fit (issue #3699)
            const Box rbx = phi[mfi].box() & geom.growPeriodicDomain(phi.nGrowVect());
            const auto rlo = lbound(rbx);
            const auto rhi = ubound(rbx);

            ParallelFor(xbx, ybx, zbx,
            // x-face
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (apx(i,j,k) == zero) {
                    if (!u_cflag(i,j,k).isCovered()) {
                        if (cflag(i,j,k).isCovered() && !cflag(i-1,j,k).isCovered()) {
                            // Extrapolate from the uncovered side, at i-1 and below
                            if (i-3 >= rlo.x) {
                                fx(i,j,k) = dxInv[0] * (p_arr(i-3,j,k) - three*p_arr(i-2,j,k) + two*p_arr(i-1,j,k));
                            } else if (i-2 >= rlo.x) {
                                fx(i,j,k) = dxInv[0] * (p_arr(i-1,j,k) - p_arr(i-2,j,k));
                            }
                        } else if (cflag(i-1,j,k).isCovered() && !cflag(i,j,k).isCovered()) {
                            // Extrapolate from the uncovered side, at i and above
                            if (i+2 <= rhi.x) {
                                fx(i,j,k) = dxInv[0] * (three*p_arr(i+1,j,k) - p_arr(i+2,j,k) - two*p_arr(i,j,k));
                            } else if (i+1 <= rhi.x) {
                                fx(i,j,k) = dxInv[0] * (p_arr(i+1,j,k) - p_arr(i,j,k));
                            }
                        }
                    }
                }
            },
            // y-face
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (apy(i,j,k) == zero) {
                    if (!v_cflag(i,j,k).isCovered()) {
                        if (cflag(i,j,k).isCovered() && !cflag(i,j-1,k).isCovered()) {
                            // Extrapolate from the uncovered side, at j-1 and below
                            if (j-3 >= rlo.y) {
                                fy(i,j,k) = dxInv[1] * (p_arr(i,j-3,k) - three*p_arr(i,j-2,k) + two*p_arr(i,j-1,k));
                            } else if (j-2 >= rlo.y) {
                                fy(i,j,k) = dxInv[1] * (p_arr(i,j-1,k) - p_arr(i,j-2,k));
                            }
                        } else if (cflag(i,j-1,k).isCovered() && !cflag(i,j,k).isCovered()) {
                            // Extrapolate from the uncovered side, at j and above
                            if (j+2 <= rhi.y) {
                                fy(i,j,k) = dxInv[1] * (three*p_arr(i,j+1,k) - p_arr(i,j+2,k) - two*p_arr(i,j,k));
                            } else if (j+1 <= rhi.y) {
                                fy(i,j,k) = dxInv[1] * (p_arr(i,j+1,k) - p_arr(i,j,k));
                            }
                        }
                    }
                }
            },
            // z-face
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (apz(i,j,k) == zero) {
                    if (!w_cflag(i,j,k).isCovered()) {
                        if (cflag(i,j,k).isCovered() && !cflag(i,j,k-1).isCovered()) {
                            // Extrapolate from the uncovered side, at k-1 and below
                            if (k-3 >= rlo.z) {
                                fz(i,j,k) = dxInv[2] * (p_arr(i,j,k-3) - three*p_arr(i,j,k-2) + two*p_arr(i,j,k-1));
                            } else if (k-2 >= rlo.z) {
                                fz(i,j,k) = dxInv[2] * (p_arr(i,j,k-1) - p_arr(i,j,k-2));
                            }
                        } else if (cflag(i,j,k-1).isCovered() && !cflag(i,j,k).isCovered()) {
                            // Extrapolate from the uncovered side, at k and above
                            if (k+2 <= rhi.z) {
                                fz(i,j,k) = dxInv[2] * (three*p_arr(i,j,k+1) - p_arr(i,j,k+2) - two*p_arr(i,j,k));
                            } else if (k+1 <= rhi.z) {
                                fz(i,j,k) = dxInv[2] * (p_arr(i,j,k+1) - p_arr(i,j,k));
                            }
                        }
                    }
                }
            });
        } // single-valued
    } // mfi
}

// Explicit template instantiations for the types we use
template void FillZeroAreaFaceFluxes(amrex::MultiFab&, amrex::Array<amrex::MultiFab,AMREX_SPACEDIM>&,
                                     const amrex::Geometry&, amrex::EBFArrayBoxFactory const&,
                                     amrex::EBFArrayBoxFactory const&,
                                     amrex::EBFArrayBoxFactory const&,
                                     amrex::EBFArrayBoxFactory const&);
template void FillZeroAreaFaceFluxes(amrex::MultiFab&, amrex::Array<amrex::MultiFab,AMREX_SPACEDIM>&,
                                     const amrex::Geometry&, amrex::EBFArrayBoxFactory const&,
                                     eb_aux_ const&,
                                     eb_aux_ const&,
                                     eb_aux_ const&);
