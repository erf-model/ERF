#include <ERF_FillPatcher.H>
#include <AMReX_Arena.H>
#include <AMReX_InterpFaceReg_3D_C.H>

using namespace amrex;

/*
 * Fill valid and ghost data with the "state data" at the given time
 *
 * @param[in] fba    BoxArray of data to be filled at fine level
 * @param[in] fdm    DistributionMapping of data to be filled at fine level
 * @param[in] fgeom  container of geometry infomation at fine level
 * @param[in] cba    BoxArray of data to be filled at coarse level
 * @param[in] cdm    DistributionMapping of data to be filled at coarse level
 * @param[in] cgeom  container of geometry infomation at coarse level
 * @param[in] nghost number of ghost cells to be filled
 * @param[in] ncomp  number of components to be filled
 * @param[in] interp interpolation operator to be used
 */

ERFFillPatcher::ERFFillPatcher (BoxArray const& fba, DistributionMapping const& fdm,
                                Geometry const& fgeom,
                                BoxArray const& cba, DistributionMapping const& cdm,
                                Geometry const& cgeom,
                                int nghost, int nghost_set,
                                int ncomp, InterpBase* interp)
    : m_fba(fba),
      m_cba(cba),
      m_fdm(fdm),
      m_cdm(cdm),
      m_fgeom(fgeom),
      m_cgeom(cgeom)
{
    AMREX_ALWAYS_ASSERT(fba.ixType() == cba.ixType());

    // Vector to hold times for coarse data
    m_crse_times.resize(2);

    // Define the coarse and fine MFs
    Define(fba, fdm, fgeom, cba, cdm, cgeom,nghost, nghost_set, ncomp, interp);
}


/*
 * Redefine the coarse and fine patch MultiFabs.
 *
 * @param[in] fba    BoxArray of data to be filled at fine level
 * @param[in] fdm    DistributionMapping of data to be filled at fine level
 * @param[in] fgeom  container of geometry infomation at fine level
 * @param[in] cba    BoxArray of data to be filled at coarse level
 * @param[in] cdm    DistributionMapping of data to be filled at coarse level
 * @param[in] cgeom  container of geometry infomation at coarse level
 * @param[in] nghost number of ghost cells to be filled
 * @param[in] ncomp  number of components to be filled
 * @param[in] interp interpolation operator to be used
 */
void ERFFillPatcher::Define (BoxArray const& fba, DistributionMapping const& fdm,
                             Geometry const& fgeom,
                             BoxArray const& cba, DistributionMapping const& cdm,
                             Geometry const& cgeom,
                             int nghost, int nghost_set,
                             int ncomp, InterpBase* interp)
{
    AMREX_ALWAYS_ASSERT(nghost <= 0);
    AMREX_ALWAYS_ASSERT(nghost_set <= 0);
    AMREX_ALWAYS_ASSERT(nghost <= nghost_set);

    // Set data memebers
    m_fba = fba; m_cba = cba;
    m_fdm = fdm; m_cdm = cdm;
    m_fgeom  = fgeom;  m_cgeom = cgeom;
    m_nghost = nghost; m_nghost_subset = nghost_set;
    m_ncomp  = ncomp;  m_interp = interp;

    // Delete old MFs if they exist
    if (m_cf_crse_data_old) m_cf_crse_data_old.reset();
    if (m_cf_crse_data_new) m_cf_crse_data_new.reset();
    if (m_cf_mask) m_cf_mask.reset();

    // Index type for the BL/BA
    IndexType m_ixt = fba.ixType();

    // Refinement ratios
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        m_ratio[idim] = m_fgeom.Domain().length(idim) / m_cgeom.Domain().length(idim);
    }

    // Coarse box list
    // NOTE: if we use face_cons_linear_interp then CoarseBox returns the grown box
    //       so we don't need to manually grow it here
    BoxList cbl;
    cbl.set(m_ixt);
    cbl.reserve(fba.size());
    for (int i(0); i < fba.size(); ++i) {
        Box coarse_box(interp->CoarseBox(fba[i], m_ratio));
        cbl.push_back(coarse_box);
    }

    // Box arrays for the coarse data
    BoxArray cf_cba(std::move(cbl));

    // Two coarse patches to hold the data to be interpolated
    m_cf_crse_data_old = std::make_unique<MultiFab> (cf_cba, fdm, m_ncomp, 0);
    m_cf_crse_data_new = std::make_unique<MultiFab> (cf_cba, fdm, m_ncomp, 0);

    // Integer masking array
    m_cf_mask = std::make_unique<iMultiFab> (fba, fdm, 1, 0);

    // Populate mask array
    if (nghost_set <= 0) {
        m_cf_mask->setVal(m_set_mask);
        BuildMask(fba,nghost_set,m_set_mask-1);
    } else {
        m_cf_mask->setVal(m_relax_mask);
        BuildMask(fba,nghost,m_relax_mask-1);
    }
}

void ERFFillPatcher::BuildMask (BoxArray const& fba,
                                int nghost,
                                int mask_val)
{
    // Minimal bounding box of fine BA plus a halo cell
    Box fba_bnd = grow(fba.minimalBox(), IntVect(1,1,1));

    // BoxList and BoxArray to store complement
    BoxList com_bl; BoxArray com_ba;

    // Compute the complement
    fba.complementIn(com_bl,fba_bnd);

    // com_bl cannot be null since we grew with halo cells
    AMREX_ALWAYS_ASSERT(com_bl.size() > 0);

    IntVect box_grow_vect(-nghost,-nghost,0);

    // cf_set_width = cf_width = 0 is a special case
    // In this case we set only the normal velocities
    // (not any cell-centered quantities) and only
    // on the coarse-fine boundary itself
    if (nghost == 0) {
        if (fba.ixType()[0] == IndexType::NODE) {
            box_grow_vect = IntVect(1,0,0);
        } else if (fba.ixType()[1] == IndexType::NODE) {
            box_grow_vect = IntVect(0,1,0);
        } else if (fba.ixType()[2] == IndexType::NODE) {
            box_grow_vect = IntVect(0,0,1);
        }
    }

    // Grow the complement boxes and trim with the bounding box
    Vector<Box>& com_bl_v = com_bl.data();
    for (int i(0); i<com_bl.size(); ++i) {
        Box& bx = com_bl_v[i];
        bx.grow(box_grow_vect);
        bx &= fba_bnd;
    }


    // Do second complement with the grown boxes
    com_ba.define(std::move(com_bl));
    com_ba.complementIn(com_bl, fba_bnd);

    // Fill mask based upon the com_bl BoxList
    for (MFIter mfi(*m_cf_mask); mfi.isValid(); ++mfi) {
        const Box& vbx = mfi.validbox();
        const Array4<int>& mask_arr = m_cf_mask->array(mfi);

        for (auto const& b : com_bl) {
            Box com_bx = vbx & b;
            ParallelFor(com_bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                mask_arr(i,j,k) = mask_val;
            });
        }
    }
}

/*
 * Register the coarse data to be used by the ERFFillPatcher
 *
 * @param[in] crse_data data at old and new time at coarse level
 * @param[in] crse_time times at which crse_data is defined
 */

void ERFFillPatcher::RegisterCoarseData (Vector<MultiFab const*> const& crse_data,
                                         Vector<Real> const& crse_time)
{
    AMREX_ALWAYS_ASSERT(crse_data.size() == 2); // old and new
    AMREX_ALWAYS_ASSERT(crse_time[1] >= crse_time[0]);

    // NOTE: CoarseBox with CellConsLinear interpolation grows the
    //       box by 1 in all directions. This pushes the domain for
    //       m_cf_crse_data into ghost cells in the z-dir. So we need
    //       to include ghost cells for crse_data when doing the copy
    IntVect src_ng = crse_data[0]->nGrowVect();
    IntVect dst_ng = m_cf_crse_data_old->nGrowVect();

    m_cf_crse_data_old->ParallelCopy(*(crse_data[0]), 0, 0, m_ncomp,
                                    src_ng, dst_ng, m_cgeom.periodicity()); // old data
    m_cf_crse_data_new->ParallelCopy(*(crse_data[1]), 0, 0, m_ncomp,
                                    src_ng, dst_ng, m_cgeom.periodicity()); // new data

    m_crse_times[0] = crse_time[0]; // time of "old" coarse data
    m_crse_times[1] = crse_time[1]; // time of "new" coarse data

    m_dt_crse = crse_time[1] - crse_time[0];
}

void ERFFillPatcher::InterpFace (MultiFab& fine,
                                 MultiFab const& crse,
                                 int mask_val)
{
    int ncomp = 1;
    IntVect ratio = m_ratio;

    FArrayBox slope;

    //
    // This box is only used to make sure we don't look outside
    //     the domain for computing the slopes in the interpolation
    // We need it to be of the type of the faces being filled
    //
    // IndexType ixt = fine.boxArray().ixType();
    // Box const& domface = convert(m_cgeom.Domain(), ixt);

    // We don't need to worry about face-based domain because this is only used in the tangential interpolation
    Box per_grown_domain = m_cgeom.Domain();
    for (int dim = 0; dim < AMREX_SPACEDIM; dim++) {
        if (m_cgeom.isPeriodic(dim)) {
            per_grown_domain.grow(dim,1);
        }
    }

    for (MFIter mfi(fine); mfi.isValid(); ++mfi)
    {
        Box const& fbx = mfi.validbox();

        slope.resize(fbx,ncomp,The_Async_Arena());

        Array4<Real> const&       fine_arr = fine.array(mfi);
        Array4<Real> const&      slope_arr = slope.array();
        Array4<Real const> const& crse_arr = crse.const_array(mfi);
        Array4<int const> const&  mask_arr = m_cf_mask->const_array(mfi);

        if (fbx.type(0) == IndexType::NODE) // x-faces
        {
            // Here do interpolation in the tangential directions
            AMREX_HOST_DEVICE_PARALLEL_FOR_3D_FLAG(RunOn::Gpu,fbx,i,j,k,
            {
                if (mask_arr(i,j,k) == mask_val) { // x-faces
                    const int ii = coarsen(i,ratio[0]);
                    if (i-ii*ratio[0] == 0) {
                        interp_face_reg(i,j,k,ratio,fine_arr,0,crse_arr,slope_arr,ncomp,per_grown_domain,0);
                    }
                }
            });

            // Here do interpolation in the normal direction
            //    using the fine values that have already been filled
            AMREX_HOST_DEVICE_PARALLEL_FOR_3D_FLAG(RunOn::Gpu,fbx,i,j,k,
            {
                if (mask_arr(i,j,k) == mask_val) {
                    const int ii = coarsen(i,ratio[0]);
                    if (i-ii*ratio[0] != 0) {
                        Real const w = static_cast<Real>(i-ii*ratio[0]) * (Real(1.)/Real(ratio[0]));
                        fine_arr(i,j,k,0) = (Real(1.)-w) * fine_arr(ii*ratio[0],j,k,0) + w * fine_arr((ii+1)*ratio[0],j,k,0);
                    }
                }
            });

        }
        else if (fbx.type(1) == IndexType::NODE) // y-faces
        {
            // Here do interpolation in the tangential directions
            AMREX_HOST_DEVICE_PARALLEL_FOR_3D_FLAG(RunOn::Gpu,fbx,i,j,k,
            {
                if (mask_arr(i,j,k) == mask_val) {
                    const int jj = coarsen(j,ratio[1]);
                    if (j-jj*ratio[1] == 0) {
                        interp_face_reg(i,j,k,ratio,fine_arr,0,crse_arr,slope_arr,ncomp,per_grown_domain,1);
                    }
                }
            });

            // Here do interpolation in the normal direction
            //    using the fine values that have already been filled
            AMREX_HOST_DEVICE_PARALLEL_FOR_3D_FLAG(RunOn::Gpu,fbx,i,j,k,
            {
                if (mask_arr(i,j,k) == mask_val) {
                    const int jj = coarsen(j,ratio[1]);
                    if (j-jj*ratio[1] != 0) {
                        Real const w = static_cast<Real>(j-jj*ratio[1]) * (Real(1.)/Real(ratio[1]));
                        fine_arr(i,j,k,0) = (Real(1.)-w) * fine_arr(i,jj*ratio[1],k,0) + w * fine_arr(i,(jj+1)*ratio[1],k,0);
                    }
                }
            });
        }
        else // z-faces
        {
            // Here do interpolation in the tangential directions
            AMREX_HOST_DEVICE_PARALLEL_FOR_3D_FLAG(RunOn::Gpu,fbx,i,j,k,
            {
                if (mask_arr(i,j,k) == mask_val) {
                    const int kk = coarsen(k,ratio[2]);
                    if (k-kk*ratio[2] == 0) {
                        interp_face_reg(i,j,k,ratio,fine_arr,0,crse_arr,slope_arr,1,per_grown_domain,2);
                    }
                }
            });

            // Here do interpolation in the normal direction
            //    using the fine values that have already been filled
            AMREX_HOST_DEVICE_PARALLEL_FOR_3D_FLAG(RunOn::Gpu,fbx,i,j,k,
            {
                if (mask_arr(i,j,k) == mask_val) {
                    const int kk = coarsen(k,ratio[2]);
                    if (k-kk*ratio[2] != 0) {
                        Real const w = static_cast<Real>(k-kk*ratio[2]) * (Real(1.)/Real(ratio[2]));
                        fine_arr(i,j,k,0) = (Real(1.)-w) * fine_arr(i,j,kk*ratio[2],0) + w * fine_arr(i,j,(kk+1)*ratio[2],0);
                    }
                }
            });
        } // IndexType::NODE
    } // MFiter
}

void ERFFillPatcher::InterpCell (MultiFab& fine,
                                 MultiFab const& crse,
                                 Vector<BCRec> const& bcr,
                                 int mask_val)
{
    int ncomp = m_ncomp;
    IntVect ratio = m_ratio;
    IndexType m_ixt = fine.boxArray().ixType();
    Box const& cdomain = convert(m_cgeom.Domain(), m_ixt);

    for (MFIter mfi(fine); mfi.isValid(); ++mfi) {
        Box const& fbx = mfi.validbox();

        Array4<Real> const&       fine_arr = fine.array(mfi);
        Array4<Real const> const& crse_arr = crse.const_array(mfi);
        Array4<int const> const&  mask_arr = m_cf_mask->const_array(mfi);

        bool run_on_gpu = Gpu::inLaunchRegion();
        amrex::ignore_unused(run_on_gpu);

        amrex::ignore_unused(m_fgeom);

        const Box& crse_region = m_interp->CoarseBox(fbx,ratio);
        Box cslope_bx(crse_region);
        for (int dim = 0; dim < AMREX_SPACEDIM; dim++) {
            if (ratio[dim] > 1) {
                cslope_bx.grow(dim,-1);
            }
        }

        FArrayBox ccfab(cslope_bx, ncomp*AMREX_SPACEDIM, The_Async_Arena());
        Array4<Real> const& tmp = ccfab.array();
        Array4<Real const> const& ctmp = ccfab.const_array();

#ifdef AMREX_USE_GPU
        AsyncArray<BCRec> async_bcr(bcr.data(), (run_on_gpu) ? ncomp : 0);
        BCRec const* bcrp = (run_on_gpu) ? async_bcr.data() : bcr.data();
#else
        BCRec const* bcrp = bcr.data();
#endif

        AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(RunOn::Gpu, cslope_bx, ncomp, i, j, k, n,
        {
            mf_cell_cons_lin_interp_mcslope(i,j,k,n, tmp, crse_arr, 0, ncomp,
                                            cdomain, ratio, bcrp);
        });

        AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(RunOn::Gpu, fbx, ncomp, i, j, k, n,
        {
            if (mask_arr(i,j,k) == mask_val) mf_cell_cons_lin_interp(i,j,k,n, fine_arr, 0, ctmp,
                                                                     crse_arr, 0, ncomp, ratio);
        });
    } // MFIter
}
