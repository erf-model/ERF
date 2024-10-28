#include <ERF.H>
#include <ERF_PhysBCFunct.H>
#include <ERF_IndexDefines.H>
#include <ERF_TimeInterpolatedData.H>
#include <ERF_FillPatcher.H>
#include <ERF_Utils.H>

using namespace amrex;

void
ERF::FillBdyCCVels (Vector<MultiFab>& mf_cc_vel)
{
    // Impose bc's at domain boundaries
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        Box domain(Geom(lev).Domain());

        int ihi = domain.bigEnd(0);
        int jhi = domain.bigEnd(1);
        int khi = domain.bigEnd(2);

        // Impose periodicity first
        mf_cc_vel[lev].FillBoundary(geom[lev].periodicity());

        for (MFIter mfi(mf_cc_vel[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            // Note that we don't fill corners here -- only the cells that share a face
            //      with interior cells -- this is all that is needed to calculate vorticity
            const Box& bx = mfi.tilebox();
            const Array4<Real>& vel_arr = mf_cc_vel[lev].array(mfi);

            if (!Geom(lev).isPeriodic(0)) {
                // Low-x side
                if (bx.smallEnd(0) <= domain.smallEnd(0)) {
                    Real mult = (phys_bc_type[0] == ERF_BC::no_slip_wall) ? -1. : 1.;
                    ParallelFor(makeSlab(bx,0,0), [=] AMREX_GPU_DEVICE(int , int j, int k) noexcept
                    {
                        vel_arr(-1,j,k,1) = mult*vel_arr(0,j,k,1); // v
                        vel_arr(-1,j,k,2) = mult*vel_arr(0,j,k,2); // w
                    });
                }

                // High-x side
                if (bx.bigEnd(0) >= domain.bigEnd(0)) {
                    Real mult = (phys_bc_type[3] == ERF_BC::no_slip_wall) ? -1. : 1.;
                    ParallelFor(makeSlab(bx,0,0), [=] AMREX_GPU_DEVICE(int , int j, int k) noexcept
                    {
                        vel_arr(ihi+1,j,k,1) = mult*vel_arr(ihi,j,k,1); // v
                        vel_arr(ihi+1,j,k,2) = mult*vel_arr(ihi,j,k,2); // w
                    });
                }
            } // !periodic

            if (!Geom(lev).isPeriodic(1)) {
                // Low-y side
                if (bx.smallEnd(1) <= domain.smallEnd(1)) {
                    Real mult = (phys_bc_type[1] == ERF_BC::no_slip_wall) ? -1. : 1.;
                    ParallelFor(makeSlab(bx,1,0), [=] AMREX_GPU_DEVICE(int i, int  , int k) noexcept
                    {
                        vel_arr(i,-1,k,0) = mult*vel_arr(i,0,k,0); // u
                        vel_arr(i,-1,k,2) = mult*vel_arr(i,0,k,2); // w
                    });
                }

                // High-y side
                if (bx.bigEnd(1) >= domain.bigEnd(1)) {
                    Real mult = (phys_bc_type[4] == ERF_BC::no_slip_wall) ? -1. : 1.;
                    ParallelFor(makeSlab(bx,1,0), [=] AMREX_GPU_DEVICE(int i, int , int k) noexcept
                    {
                        vel_arr(i,jhi+1,k,0) = mult*vel_arr(i,jhi,k,0); // u
                        vel_arr(i,jhi+1,k,2) = mult*-vel_arr(i,jhi,k,2); // w
                    });
                }
            } // !periodic

            if (!Geom(lev).isPeriodic(2)) {
                // Low-z side
                if (bx.smallEnd(2) <= domain.smallEnd(2)) {
                    Real mult = (phys_bc_type[2] == ERF_BC::no_slip_wall) ? -1. : 1.;
                    ParallelFor(makeSlab(bx,2,0), [=] AMREX_GPU_DEVICE(int i, int j, int) noexcept
                    {
                        vel_arr(i,j,-1,0) = mult*vel_arr(i,j,0,0); // u
                        vel_arr(i,j,-1,1) = mult*vel_arr(i,j,0,1); // v
                    });
                }

                // High-z side
                if (bx.bigEnd(2) >= domain.bigEnd(2)) {
                    Real mult = (phys_bc_type[5] == ERF_BC::no_slip_wall) ? -1. : 1.;
                    ParallelFor(makeSlab(bx,2,0), [=] AMREX_GPU_DEVICE(int i, int j, int) noexcept
                    {
                        vel_arr(i,j,khi+1,0) = mult*vel_arr(i,j,khi,0); // u
                        vel_arr(i,j,khi+1,1) = mult*vel_arr(i,j,khi,1); // v
                    });
                }
            } // !periodic
        } // MFIter

    } // lev
}
