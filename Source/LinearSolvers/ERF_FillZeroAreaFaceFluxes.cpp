#include "ERF.H"
#include "ERF_Utils.H"

using namespace amrex;

/**
 * Compute phi gradients where the area of the face is zero
 */
void
FillZeroAreaFaceFluxes (MultiFab& phi, Array<MultiFab,AMREX_SPACEDIM>& fluxes,
                        const Geometry& geom, EBFArrayBoxFactory const& ebfact,
                        eb_aux_ const& ebfact_u, eb_aux_ const& ebfact_v, eb_aux_ const& ebfact_w)
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

            ParallelFor(xbx, ybx, zbx,
            // x-face
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (apx(i,j,k) == Real(0.0)) {
                    if (!u_cflag(i,j,k).isCovered()) {
                        if (cflag(i,j,k).isCovered() && !cflag(i-1,j,k).isCovered()) {
                            fx(i,j,k) = dxInv[0] * (p_arr(i-3,j,k) - 3.*p_arr(i-2,j,k) + 2.*p_arr(i-1,j,k));
                        } else if (cflag(i-1,j,k).isCovered() && !cflag(i,j,k).isCovered()) {
                            fx(i,j,k) = dxInv[0] * (3.*p_arr(i+1,j,k) - p_arr(i+2,j,k) - 2.*p_arr(i,j,k));
                        }
                    }
                }
            },
            // y-face
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (apy(i,j,k) == Real(0.0)) {
                    if (!v_cflag(i,j,k).isCovered()) {
                        if (cflag(i,j,k).isCovered() && !cflag(i,j-1,k).isCovered()) {
                            fy(i,j,k) = dxInv[1] * (p_arr(i,j-3,k) - 3.*p_arr(i,j-2,k) + 2.*p_arr(i,j-1,k));
                        } else if (cflag(i,j-1,k).isCovered() && !cflag(i,j,k).isCovered()) {
                            fy(i,j,k) = dxInv[1] * (3.*p_arr(i,j+1,k) - p_arr(i,j+2,k) - 2.*p_arr(i,j,k));
                        }
                    }
                }
            },
            // z-face
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (apz(i,j,k) == Real(0.0)) {
                    if (!w_cflag(i,j,k).isCovered()) {
                        if (cflag(i,j,k).isCovered() && !cflag(i,j,k-1).isCovered()) {
                            fz(i,j,k) = dxInv[2] * (p_arr(i,j,k-3) - 3.*p_arr(i,j,k-2) + 2.*p_arr(i,j,k-1));
                        } else if (cflag(i,j,k-1).isCovered() && !cflag(i,j,k).isCovered()) {
                            fz(i,j,k) = dxInv[2] * (3.*p_arr(i,j,k+1) - p_arr(i,j,k+2) - 2.*p_arr(i,j,k));
                        }
                    }
                }
            });
        } // single-valued
    } // mfi
}
