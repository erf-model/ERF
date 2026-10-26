#include <ERF.H>
#include <ERF_PhysBCFunct.H>
#include <ERF_IndexDefines.H>
#include <ERF_TimeInterpolatedData.H>
#include <ERF_FillPatcher.H>
#include <ERF_Utils.H>

using namespace amrex;

/**
 * Fill boundary ghost cells for cell-centered velocities.
 *
 * @param[in,out] mf_cc_vel cell-centered velocity MultiFab to be filled
 * @param[in]     lev_geom  geometry for the level being filled
 */
void
ERF::FillBdyCCVels (MultiFab& mf_cc_vel, Geometry& lev_geom)
{
    // Impose bc's at domain boundaries
    Box domain(lev_geom.Domain());

    // Each side fills the one ghost layer outside the domain from the domain's edge cell
    int ilo = domain.smallEnd(0), ihi = domain.bigEnd(0);
    int jlo = domain.smallEnd(1), jhi = domain.bigEnd(1);
    int klo = domain.smallEnd(2), khi = domain.bigEnd(2);

    // Impose periodicity first
    mf_cc_vel.FillBoundary(lev_geom.periodicity());

    int jper = (lev_geom.isPeriodic(1));
    int kper = (lev_geom.isPeriodic(2));

    for (MFIter mfi(mf_cc_vel, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        const Array4<Real>& vel_arr = mf_cc_vel.array(mfi);

        if (!lev_geom.isPeriodic(0)) {
            // Low-x side
            if (bx.smallEnd(0) <= ilo) {
                Real multn = ( (phys_bc_type[0] == ERF_BC::slip_wall   ) ||
                               (phys_bc_type[0] == ERF_BC::no_slip_wall) ||
                               (phys_bc_type[0] == ERF_BC::symmetry    ) ) ? -one : one;
                Real multt = (phys_bc_type[0] == ERF_BC::no_slip_wall) ? -one : one;
                Box gbx(bx); gbx.grow(1,jper); gbx.grow(2,kper);
                ParallelFor(makeSlab(gbx,0,ilo-1), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                {
                    vel_arr(i,j,k,0) = multn*vel_arr(ilo,j,k,0); // u
                    vel_arr(i,j,k,1) = multt*vel_arr(ilo,j,k,1); // v
                    vel_arr(i,j,k,2) = multt*vel_arr(ilo,j,k,2); // w
                });
            }

            // High-x side
            if (bx.bigEnd(0) >= domain.bigEnd(0)) {
                Real multn = ( (phys_bc_type[3] == ERF_BC::slip_wall   ) ||
                               (phys_bc_type[3] == ERF_BC::no_slip_wall) ||
                               (phys_bc_type[3] == ERF_BC::symmetry    ) ) ? -one : one;
                Real multt = (phys_bc_type[3] == ERF_BC::no_slip_wall) ? -one : one;
                Box gbx(bx); gbx.grow(1,jper); gbx.grow(2,kper);
                ParallelFor(makeSlab(gbx,0,ihi+1), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                {
                    vel_arr(i,j,k,0) = multn*vel_arr(ihi,j,k,0); // u
                    vel_arr(i,j,k,1) = multt*vel_arr(ihi,j,k,1); // v
                    vel_arr(i,j,k,2) = multt*vel_arr(ihi,j,k,2); // w
                });
            }
        } // !periodic

        if (!lev_geom.isPeriodic(1)) {
            // Low-y side
            if (bx.smallEnd(1) <= domain.smallEnd(1)) {
                Real multn = ( (phys_bc_type[1] == ERF_BC::slip_wall   ) ||
                               (phys_bc_type[1] == ERF_BC::no_slip_wall) ||
                               (phys_bc_type[1] == ERF_BC::symmetry    ) ) ? -one : one;
                Real multt = (phys_bc_type[1] == ERF_BC::no_slip_wall) ? -one : one;
                Box gbx(bx); gbx.grow(0,1); gbx.grow(2,kper);
                ParallelFor(makeSlab(gbx,1,jlo-1), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                {
                    vel_arr(i,j,k,0) = multt*vel_arr(i,jlo,k,0); // u
                    vel_arr(i,j,k,1) = multn*vel_arr(i,jlo,k,1); // v
                    vel_arr(i,j,k,2) = multt*vel_arr(i,jlo,k,2); // w
                });
            }

            // High-y side
            if (bx.bigEnd(1) >= domain.bigEnd(1)) {
                Real multn = ( (phys_bc_type[4] == ERF_BC::slip_wall   ) ||
                               (phys_bc_type[4] == ERF_BC::no_slip_wall) ||
                               (phys_bc_type[4] == ERF_BC::symmetry    ) ) ? -one : one;
                Real multt = (phys_bc_type[4] == ERF_BC::no_slip_wall) ? -one : one;
                Box gbx(bx); gbx.grow(0,1); gbx.grow(2,kper);
                ParallelFor(makeSlab(gbx,1,jhi+1), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                {
                    vel_arr(i,j,k,0) = multt*vel_arr(i,jhi,k,0); // u
                    vel_arr(i,j,k,1) = multn*vel_arr(i,jhi,k,1); // v
                    vel_arr(i,j,k,2) = multt*vel_arr(i,jhi,k,2); // w
                });
            }
        } // !periodic

        if (!lev_geom.isPeriodic(2)) {
            // Low-z side
            if (bx.smallEnd(2) <= domain.smallEnd(2)) {
                Real multn = ( (phys_bc_type[2] == ERF_BC::slip_wall   ) ||
                               (phys_bc_type[2] == ERF_BC::no_slip_wall) ||
                               (phys_bc_type[2] == ERF_BC::symmetry    ) ) ? -one : one;
                Real multt = (phys_bc_type[2] == ERF_BC::no_slip_wall) ? -one : one;
                Box gbx(bx); gbx.grow(0,1); gbx.grow(1,1);
                ParallelFor(makeSlab(gbx,2,klo-1), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                {
                    vel_arr(i,j,k,0) = multt*vel_arr(i,j,klo,0); // u
                    vel_arr(i,j,k,1) = multt*vel_arr(i,j,klo,1); // v
                    vel_arr(i,j,k,2) = multn*vel_arr(i,j,klo,2); // w
                });
            }

            // High-z side
            if (bx.bigEnd(2) >= domain.bigEnd(2)) {
                Real multn = ( (phys_bc_type[5] == ERF_BC::slip_wall   ) ||
                               (phys_bc_type[5] == ERF_BC::no_slip_wall) ||
                               (phys_bc_type[5] == ERF_BC::symmetry    ) ) ? -one : one;
                Real multt = (phys_bc_type[5] == ERF_BC::no_slip_wall) ? -one : one;
                Box gbx(bx); gbx.grow(0,1); gbx.grow(1,1);
                ParallelFor(makeSlab(gbx,2,khi+1), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                {
                    vel_arr(i,j,k,0) = multt*vel_arr(i,j,khi,0); // u
                    vel_arr(i,j,k,1) = multt*vel_arr(i,j,khi,1); // v
                    vel_arr(i,j,k,2) = multn*vel_arr(i,j,khi,2); // w
                });
            }
        } // !periodic
    } // MFIter

    // Impose periodicity again
    mf_cc_vel.FillBoundary(lev_geom.periodicity());
}
