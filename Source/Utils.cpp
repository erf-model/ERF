/** \addtogroup Utilities
 * @{
 *
 * \file Utils.cpp
 *
 */

#include <AMReX_MultiFabUtil.H>
#include <AMReX_BCRec.H>
#include <utils.H>
#include <utils_K.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_FillPatchUtil.H>

using namespace amrex;

void
create_umac_grown (int lev, int nGrow, BoxArray& fine_grids,
                   const Geometry& crse_geom,
                   const Geometry& fine_geom,
                   const Array<MultiFab*,AMREX_SPACEDIM> u_mac_crse,
                   const Array<MultiFab*,AMREX_SPACEDIM> u_mac_fine,
                   const IntVect& crse_ratio)
{
    BL_PROFILE("create_umac_grown()");

    if (lev > 0)
    {
        BoxList bl = amrex::GetBndryCells(fine_grids,nGrow);

        BoxArray f_bnd_ba(std::move(bl));

        BoxArray c_bnd_ba = f_bnd_ba; c_bnd_ba.coarsen(crse_ratio);

        c_bnd_ba.maxSize(32);

        f_bnd_ba = c_bnd_ba; f_bnd_ba.refine(crse_ratio);

        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        {
            // We must make sure the coarse data has the periodic
            //    boundaries filled or we will be operating below
            //    on uninitialized data
            u_mac_crse[idim]->FillBoundary(crse_geom.periodicity());

            //
            // crse_src & fine_src must have same parallel distribution.
            // We'll use the KnapSack distribution for the fine_src_ba.
            // Since fine_src_ba should contain more points, this'll lead
            // to a better distribution.
            //
            BoxArray crse_src_ba(c_bnd_ba), fine_src_ba(f_bnd_ba);

            crse_src_ba.surroundingNodes(idim);
            fine_src_ba.surroundingNodes(idim);

            const int N = fine_src_ba.size();

            std::vector<amrex::Long> wgts(N);

#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (int i = 0; i < N; i++)
                wgts[i] = fine_src_ba[i].numPts();

            DistributionMapping dm;
            // This DM won't be put into the cache.
            dm.KnapSackProcessorMap(wgts,ParallelDescriptor::NProcs());
            /* Compiling on Windows fails due to line above. The message is:
             *
             * D:\a\ERF\ERF\Source\Utils.cpp(69,70): error C2664:
             * 'void amrex::DistributionMapping::KnapSackProcessorMap(const
             * std::vector<amrex::Long,std::allocator<amrex::Long>> &,int,amrex::Real *,bool,int,bool)':
             * cannot convert argument 1 from 'std::vector<long,std::allocator<long>>' to
             * 'const std::vector<amrex::Long,std::allocator<amrex::Long>> &'
             * [D:\a\ERF\ERF\Build\Exec\erf_srclib.vcxproj]
             */


            // FIXME
            // Declaring in this way doesn't work. I think it's because the box arrays
            // have been changed and each src box is not completely contained within a
            // single box in the Factory's BA
            // For now, coarse-fine boundary doesn't intersect EB, so should be okay...
            // MultiFab crse_src(crse_src_ba, dm, 1, 0, MFInfo(), getLevel(lev-1).Factory());
            // MultiFab fine_src(fine_src_ba, dm, 1, 0, MFInfo(), Factory());
            MultiFab crse_src(crse_src_ba, dm, 1, 0);
            MultiFab fine_src(fine_src_ba, dm, 1, 0);

            crse_src.setVal(1.e200);
            fine_src.setVal(1.e200);
            //
            // We want to fill crse_src from lower level u_mac including u_mac's grow cells.
            //
            const MultiFab& u_macLL = *u_mac_crse[idim];
            crse_src.ParallelCopy(u_macLL,0,0,1,u_macLL.nGrow(),0);

            const amrex::GpuArray<int,AMREX_SPACEDIM> c_ratio = {crse_ratio[0],crse_ratio[1],crse_ratio[2]};

            //
            // Fill fine values with piecewise-constant interp of coarse data.
            // Operate only on faces that overlap--ie, only fill the fine faces that make up each
            // coarse face, leave the in-between faces alone.
            //
#ifdef _OPENMP
#pragma omp parallel
#endif
            for (MFIter mfi(crse_src); mfi.isValid(); ++mfi)
            {
                const Box& box       = crse_src[mfi].box();
                auto const& crs_arr  = crse_src.array(mfi);
                auto const& fine_arr = fine_src.array(mfi);

                ParallelFor(box,[crs_arr,fine_arr,idim,c_ratio]
                AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                   int idx[3] = {i*c_ratio[0],j*c_ratio[1],k*c_ratio[2]};
                   // dim1 and dim2 are the complements of idim
                   int dim1 = ( idim != 0 ) ? 0 : 1 ;
                   int dim2 = ( idim != 0 ) ? ( ( idim == 2 ) ? 1 : 2 ) : 2 ;
                   for (int n1 = 0; n1 < c_ratio[dim1]; n1++) {
                      for (int n2 = 0; n2 < c_ratio[dim2]; n2++) {
                         int id[3] = {idx[0],idx[1],idx[2]};
                         id[dim1] += n1;
                         id[dim2] += n2;
                         fine_arr(id[0],id[1],id[2]) = crs_arr(i,j,k);
                      }
                   }
                });
            }
            crse_src.clear();
            //
            // Replace pc-interpd fine data with preferred u_mac data at
            // this level u_mac valid only on surrounding faces of valid
            // region - this op will not fill grow region.
            //
            fine_src.ParallelCopy(*u_mac_fine[idim]);

            //
            // Interpolate unfilled grow cells using best data from
            // surrounding faces of valid region, and pc-interpd data
            // on fine faces overlaying coarse edges.
            //
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(fine_src); mfi.isValid(); ++mfi)
            {
                const int  nComp = 1;
                const Box& fbox  = fine_src[mfi].box();
                auto const& fine_arr = fine_src.array(mfi);

                if (fbox.type(0) == IndexType::NODE)
                {
                  AMREX_HOST_DEVICE_PARALLEL_FOR_4D(fbox,nComp,i,j,k,n,
                  {
                    face_interp_x(i,j,k,n,fine_arr,c_ratio);
                  });
                }
                else if (fbox.type(1) == IndexType::NODE)
                {
                  AMREX_HOST_DEVICE_PARALLEL_FOR_4D(fbox,nComp,i,j,k,n,
                  {
                    face_interp_y(i,j,k,n,fine_arr,c_ratio);
                  });
                }
                else
                {
                  AMREX_HOST_DEVICE_PARALLEL_FOR_4D(fbox,nComp,i,j,k,n,
                          {
                    face_interp_z(i,j,k,n,fine_arr,c_ratio);
                  });
                }
            } // mfi

            MultiFab u_mac_save(u_mac_fine[idim]->boxArray(),u_mac_fine[idim]->DistributionMap(),1,0,MFInfo(),
                                u_mac_fine[idim]->Factory());
            u_mac_save.ParallelCopy(*u_mac_fine[idim]);
            u_mac_fine[idim]->ParallelCopy(fine_src,0,0,1,0,nGrow);
            u_mac_fine[idim]->ParallelCopy(u_mac_save);

        } // idim
    } // lev > 0

    for (int n = 0; n < AMREX_SPACEDIM; ++n)
    {
        u_mac_fine[n]->FillBoundary(fine_geom.periodicity());
    }
}

/** @}*/
