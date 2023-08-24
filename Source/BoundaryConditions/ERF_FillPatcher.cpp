#include <ERF_FillPatcher.H>

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

ERFFillPatcher::ERFFillPatcher (BoxArray const& fba, DistributionMapping  fdm,
                                Geometry const& fgeom,
                                BoxArray  cba, DistributionMapping  cdm,
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
    Define(fba, fdm, fgeom, cba, cdm, cgeom,
           nghost, nghost_set, ncomp, interp);
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
void ERFFillPatcher::Define (BoxArray const& fba, DistributionMapping  fdm,
                             Geometry const& fgeom,
                             BoxArray  cba, DistributionMapping  cdm,
                             Geometry const& cgeom,
                             int nghost, int nghost_set,
                             int ncomp, InterpBase* interp)
{
    AMREX_ALWAYS_ASSERT(nghost < 0);
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
    BoxList cbl;
    cbl.set(m_ixt);
    cbl.reserve(fba.size());
    for (int i(0); i<fba.size(); ++i) {
        cbl.push_back(interp->CoarseBox(fba[i], m_ratio));
    }

    // Box arrays for the coarse data
    BoxArray cf_cba(std::move(cbl));
    DistributionMapping cf_dm(cf_cba);

    // Two coarse patches to hold the data to be interpolated
    m_cf_crse_data_old = std::make_unique<MultiFab> (cf_cba, fdm, m_ncomp, 0);
    m_cf_crse_data_new = std::make_unique<MultiFab> (cf_cba, fdm, m_ncomp, 0);

    // Integer masking array
    m_cf_mask = std::make_unique<iMultiFab> (fba, fdm, 1, 0);
    m_cf_mask->setVal(0);

    // Populate mask array
    BuildMask(fba,nghost,m_relax_mask);
    if (nghost_set < 0) BuildMask(fba,nghost_set,m_set_mask);
}

void ERFFillPatcher::BuildMask (BoxArray const& fba,
                                int nghost,
                                int mask_val)
{
    // Index type for the BA
    IndexType m_ixt = fba.ixType();

    // Get interior halo cells for every box in fine ba
    Vector<BoxList> halo_v; halo_v.resize(fba.size());
    for (int ibox(0); ibox<halo_v.size(); ++ibox) {
        BoxList& halo = halo_v[ibox];
        halo.set(m_ixt);
        Box vbx = fba[ibox];
        Box sbx = amrex::grow(vbx, IntVect(nghost,nghost,0));
        BoxList const& bndry = amrex::boxDiff(vbx, sbx);
        if (bndry.isNotEmpty()) {
            halo.join(bndry);
        }
    }

    // Get interface unions (extend into boxes by nghost in normal dir)
    Vector<Vector<Box>> fba_un;
    fba_un.resize(fba.size());
    for (int ibox(0); ibox<fba_un.size()-1; ++ibox) {
        const Box& src_bx = fba[ibox];
        const IntVect& src_se = src_bx.smallEnd();
        const IntVect& src_be = src_bx.bigEnd();
        for (int jbox(ibox+1); jbox<fba_un.size(); ++jbox) {
            const Box& dst_bx = fba[jbox];
            const IntVect& dst_se = dst_bx.smallEnd();
            const IntVect& dst_be = dst_bx.bigEnd();
            for (int idim(0); idim < AMREX_SPACEDIM-1; ++idim) {
                Box src_gbx_lo = amrex::growLo(src_bx, idim, 1);
                Box src_gbx_hi = amrex::growHi(src_bx, idim, 1);
                if (src_gbx_lo.intersects(dst_bx)) {
                    IntVect src_u_se = src_se; IntVect src_u_be = src_be;
                    IntVect dst_u_se = dst_se; IntVect dst_u_be = dst_be;

                    int jdim = (idim==0) ? 1 : 0;
                    src_u_be[idim] = src_u_se[idim] - nghost - 1; // nghost is negative
                    src_u_se[jdim] = std::max(src_se[jdim], dst_se[jdim]);
                    src_u_be[jdim] = std::min(src_be[jdim], dst_be[jdim]);

                    dst_u_se[idim] = dst_u_be[idim] + nghost + 1; // nghost is negative
                    dst_u_se[jdim] = src_u_se[jdim];
                    dst_u_be[jdim] = src_u_be[jdim];

                    // Test if end points lie in any other box
                    bool src_lo_found = false; bool src_hi_found = false;
                    bool dst_lo_found = false; bool dst_hi_found = false;
                    IntVect src_u_gse = src_u_se; src_u_gse[jdim] -= 1;
                    IntVect src_u_gbe = src_u_be; src_u_gbe[jdim] += 1;
                    IntVect dst_u_gse = dst_u_se; dst_u_gse[jdim] -= 1;
                    IntVect dst_u_gbe = dst_u_be; dst_u_gbe[jdim] += 1;
                    for (int kbox(0); kbox<fba.size(); ++kbox) {
                        Box test_bx = fba[kbox];
                        if (test_bx.contains(src_u_gse)) src_lo_found = true;
                        if (test_bx.contains(dst_u_gse)) dst_lo_found = true;
                        if (test_bx.contains(src_u_gbe)) src_hi_found = true;
                        if (test_bx.contains(dst_u_gbe)) dst_hi_found = true;
                    }
                    if (!src_lo_found || !dst_lo_found) {src_u_se[jdim] -= nghost; dst_u_se[jdim] -= nghost;}
                    if (!src_hi_found || !dst_hi_found) {src_u_be[jdim] += nghost; dst_u_be[jdim] += nghost;}

                    // Build the face union box and add to BL
                    Box src_u_bx(src_u_se,src_u_be,m_ixt);
                    Box dst_u_bx(dst_u_se,dst_u_be,m_ixt);
                    fba_un[ibox].push_back(src_u_bx);
                    fba_un[jbox].push_back(dst_u_bx);

                } else if (src_gbx_hi.intersects(dst_bx)) {
                    IntVect src_u_se = src_se; IntVect src_u_be = src_be;
                    IntVect dst_u_se = dst_se; IntVect dst_u_be = dst_be;

                    int jdim = (idim==0) ? 1 : 0;
                    src_u_se[idim] = src_u_be[idim] + nghost + 1; // nghost is negative
                    src_u_se[jdim] = std::max(src_se[jdim], dst_se[jdim]);
                    src_u_be[jdim] = std::min(src_be[jdim], dst_be[jdim]);

                    dst_u_be[idim] = dst_u_se[idim] - nghost - 1; // nghost is negative
                    dst_u_se[jdim] = src_u_se[jdim];
                    dst_u_be[jdim] = src_u_be[jdim];

                    // Test if end points lie in any other box
                    bool src_lo_found = false; bool src_hi_found = false;
                    bool dst_lo_found = false; bool dst_hi_found = false;
                    IntVect src_u_gse = src_u_se; src_u_gse[jdim] -= 1;
                    IntVect src_u_gbe = src_u_be; src_u_gbe[jdim] += 1;
                    IntVect dst_u_gse = dst_u_se; dst_u_gse[jdim] -= 1;
                    IntVect dst_u_gbe = dst_u_be; dst_u_gbe[jdim] += 1;
                    for (int kbox(0); kbox<fba.size(); ++kbox) {
                        Box test_bx = fba[kbox];
                        if (test_bx.contains(src_u_gse)) src_lo_found = true;
                        if (test_bx.contains(dst_u_gse)) dst_lo_found = true;
                        if (test_bx.contains(src_u_gbe)) src_hi_found = true;
                        if (test_bx.contains(dst_u_gbe)) dst_hi_found = true;
                    }
                    if (!src_lo_found || !dst_lo_found) {src_u_se[jdim] -= nghost; dst_u_se[jdim] -= nghost;}
                    if (!src_hi_found || !dst_hi_found) {src_u_be[jdim] += nghost; dst_u_be[jdim] += nghost;}

                    // Build the face union box and add to BL
                    Box src_u_bx(src_u_se,src_u_be,m_ixt);
                    Box dst_u_bx(dst_u_se,dst_u_be,m_ixt);
                    fba_un[ibox].push_back(src_u_bx);
                    fba_un[jbox].push_back(dst_u_bx);
                } // intersects
            } // idim
        } // jbox
    } // ibox

    // Subtract off the interface regions from interior halo cells
    Vector<BoxList> mod_halo_v; mod_halo_v.resize(fba.size());
    for (int i(0); i<mod_halo_v.size(); ++i) { mod_halo_v[i].set(m_ixt); }

    for (int ibox(0); ibox<halo_v.size(); ++ibox) {
        Vector<Box>& m_halo = halo_v[ibox].data();
        BoxList&   mod_halo = mod_halo_v[ibox];
        for (int ihalo(0); ihalo<m_halo.size(); ++ihalo) {
            Box& halo_bx = m_halo[ihalo];
            BoxList tmp; tmp.set(m_ixt); tmp.push_back(halo_bx);
            Vector<Box>& m_tmp = tmp.data();
            for (int ubox(0); ubox<fba_un[ibox].size(); ++ubox) {
                Box& src_u_bx = fba_un[ibox][ubox];
                for (int dbox(0); dbox<m_tmp.size(); ++dbox) {
                    Box test = m_tmp[dbox];
                    if (test.intersects(src_u_bx)) {
                        Box halo_u_bx = test & src_u_bx;
                        BoxList const& bndry = amrex::boxDiff(test, halo_u_bx);
                        m_tmp.erase(m_tmp.begin() + dbox);
                        if (bndry.isNotEmpty()) tmp.join(bndry);
                    } // intersects
                } // dbox
            } // ubox
            mod_halo.join(tmp);
        } // ihalo
    } // ibox

    // Fill mask based upon the mod_halo BoxList
    for (MFIter mfi(*m_cf_mask); mfi.isValid(); ++mfi) {
        const Box& tbx = mfi.tilebox();
        const Array4<int>& mask_arr = m_cf_mask->array(mfi);
        BoxList& mod_halo = mod_halo_v[mfi.index()];
        const Vector<Box>& m_halo = mod_halo.data();

        for (int ibox(0); ibox<m_halo.size(); ++ibox) {
            Box mbx = tbx & m_halo[ibox];
            amrex::ParallelFor(mbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
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
                                 MultiFab& crse,
                                 int mask_val)
{
  int ncomp = m_ncomp;
  IntVect ratio = m_ratio;

  for (MFIter mfi(fine); mfi.isValid(); ++mfi) {
      Box const& fbx = mfi.validbox();

      Array4<Real> const&       fine_arr = fine.array(mfi);
      Array4<Real const> const& crse_arr = crse.const_array(mfi);
      Array4<int const> const&  mask_arr = m_cf_mask->const_array(mfi);

      if (fbx.type(0) == IndexType::NODE) {
          AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(RunOn::Gpu,fbx,ncomp,i,j,k,n,
          {
              if (mask_arr(i,j,k) == mask_val) face_linear_interp_x(i,j,k,n,fine_arr,crse_arr,ratio);
          });
      } else if (fbx.type(1) == IndexType::NODE) {
          AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(RunOn::Gpu,fbx,ncomp,i,j,k,n,
          {
              if (mask_arr(i,j,k) == mask_val) face_linear_interp_y(i,j,k,n,fine_arr,crse_arr,ratio);
          });
      } else {
          AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(RunOn::Gpu,fbx,ncomp,i,j,k,n,
          {
              if (mask_arr(i,j,k) == mask_val) face_linear_interp_z(i,j,k,n,fine_arr,crse_arr,ratio);
          });
      } // IndexType::NODE
  } // MFiter
}

void ERFFillPatcher::InterpCell (MultiFab& fine,
                                 MultiFab& crse,
                                 Vector<BCRec> const& bcr,
                                 int mask_val)
{

  int ncomp = m_ncomp;
  IntVect ratio = m_ratio;
  IndexType m_ixt = fine.boxArray().ixType();
  Box const& cdomain = amrex::convert(m_cgeom.Domain(), m_ixt);

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

      FArrayBox ccfab(cslope_bx, ncomp*AMREX_SPACEDIM);
      Array4<Real> const& tmp = ccfab.array();
      Array4<Real const> const& ctmp = ccfab.const_array();

#ifdef AMREX_USE_GPU
      AsyncArray<BCRec> async_bcr(bcr.data(), (run_on_gpu) ? ncomp : 0);
      BCRec const* bcrp = (run_on_gpu) ? async_bcr.data() : bcr.data();

      Elixir cceli;
      if (run_on_gpu) { cceli = ccfab.elixir(); }
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
