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
                                int nghost, int nghost_subset,
                                int ncomp, InterpBase* interp)
    : m_fba(fba),
      m_cba(std::move(cba)),
      m_fdm(std::move(fdm)),
      m_cdm(std::move(cdm)),
      m_fgeom(fgeom),
      m_cgeom(cgeom),
      m_nghost(nghost),
      m_nghost_subset(nghost_subset),
      m_ncomp(ncomp),
      m_interp(interp)
{
    AMREX_ALWAYS_ASSERT(nghost < 0);
    AMREX_ALWAYS_ASSERT(nghost_subset < 0);
    AMREX_ALWAYS_ASSERT(nghost < nghost_subset);
    AMREX_ALWAYS_ASSERT(fba.ixType() == cba.ixType());

    // Vector to hold times for coarse data
    m_crse_times.resize(2);

    // Index type for the BL/BA
    IndexType m_ixt = fba.ixType();

    // Box bounding all the fine boxes
    Box const fine_valid_box = amrex::grow(fba.minimalBox(), IntVect(nghost,nghost,0));
    Box const fine_valid_box_subset = amrex::grow(fba.minimalBox(), IntVect(nghost_subset,nghost_subset,0));

    // Refinement ratios
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        m_ratio[idim] = m_fgeom.Domain().length(idim) / m_cgeom.Domain().length(idim);
    }

    // Fine box list (interior halo regions exactly)
    BoxList bl;
    bl.set(m_ixt);
    for (int ibox = 0; ibox < fba.size(); ++ibox) {
        Box vbx = fba[ibox];
        BoxList const& bndry = amrex::boxDiff(vbx, fine_valid_box);
        if (bndry.isNotEmpty()) {
            bl.join(bndry);
        }
    }

    // Fine subset box list (interior halo regions exactly)
    BoxList bl_subset;
    bl_subset.set(m_ixt);
    for (int ibox = 0; ibox < fba.size(); ++ibox) {
        Box vbx = fba[ibox];
        BoxList const& bndry = amrex::boxDiff(vbx, fine_valid_box_subset);
        if (bndry.isNotEmpty()) {
            bl_subset.join(bndry);
        }
    }

    // Grown fine box list (expand interior halo to touch crse face)
    BoxList gbl;
    gbl.set(m_ixt);
    for (auto bl_iter = bl.begin(); bl_iter < bl.end(); bl_iter++) {
      const int* lo = (*bl_iter).loVect();
      const int* hi = (*bl_iter).hiVect();
      IntVect glo(0); IntVect ghi(0);
      for (int idim(0); idim < AMREX_SPACEDIM; ++idim) {
        int rlo = lo[idim]%m_ratio[idim];
        int rhi = hi[idim]%m_ratio[idim];
        glo[idim] = lo[idim]; ghi[idim] = hi[idim];
        if (rlo > 0) glo[idim] = lo[idim] - rlo;
        if (rhi > 0) ghi[idim] = hi[idim] - rhi + m_ratio[idim];
      }
      Box gbx(glo,ghi,m_ixt);
      gbl.push_back(gbx);
    }

    // Coarse box list (coincides with grown fine box list)
    BoxList cbl;
    cbl.set(m_ixt);
    cbl.reserve(gbl.size());
    for (auto const& b : gbl) {
        cbl.push_back(interp->CoarseBox(b, m_ratio));
    }

    // Box arrays for fine and coarse
    BoxArray gcf_fba(std::move(gbl));
    BoxArray gcf_cba(std::move(cbl));
    BoxArray cf_fba_s(std::move(bl_subset));
    DistributionMapping gcf_dm(gcf_fba);
    DistributionMapping cf_dm_s(cf_fba_s);

    /*
    amrex::Print() << "\n";
    amrex::Print() << "BA: " << gcf_fba << ' ' << gcf_cba << "\n";
    amrex::Print() << "DM: " << gcf_dm << "\n";
    amrex::Print() << "BA_s: " << cf_fba_s << "\n";
    amrex::Print() << "DM_s: " << cf_dm_s << "\n";
    amrex::Print() << "\n";
    */

    // Fine patch to hold the time-interpolated state
    m_cf_fine_data.define(gcf_fba, gcf_dm, m_ncomp, 0);

    // Fine subset patch to hold the time-interpolated state
    //m_cf_fine_subset_data.define(cf_fba_s, cf_dm_s, m_ncomp, 0);
    m_cf_fine_subset_data.define(cf_fba_s, gcf_dm, m_ncomp, 0);

    // Two coarse patches to hold the data to be interpolated
    m_cf_crse_data.emplace_back(gcf_cba, gcf_dm, m_ncomp, 0);
    m_cf_crse_data.emplace_back(gcf_cba, gcf_dm, m_ncomp, 0);
}


/*
 * Register the coarse data to be used by the ERFFillPatcher
 *
 * @param[in] crse_data data at old and new time at coarse level
 * @param[in] crse_time times at which crse_data is defined
 */

void ERFFillPatcher::registerCoarseData (Vector<MultiFab const*> const& crse_data,
                                         Vector<Real> const& crse_time)
{
    AMREX_ALWAYS_ASSERT(crse_data.size() == 2); // old and new
    AMREX_ALWAYS_ASSERT(crse_time[1] >= crse_time[0]);

    // NOTE: CoarseBox with CellConsLinear interpolation grows the
    //       box by 1 in all directions. This pushes the domain for
    //       m_cf_crse_data into ghost cells in the z-dir. So we need
    //       to include ghost cells for crse_data when doing the copy
    IntVect src_ng = crse_data[0]->nGrowVect();
    IntVect dst_ng = m_cf_crse_data[0].nGrowVect();
    m_cf_crse_data[0].ParallelCopy(*crse_data[0], 0, 0, m_ncomp,
                                   src_ng, dst_ng, m_cgeom.periodicity()); // old data
    m_cf_crse_data[1].ParallelCopy(*crse_data[1], 0, 0, m_ncomp,
                                   src_ng, dst_ng, m_cgeom.periodicity()); // new data

    m_crse_times[0] = crse_time[0]; // time of "old" coarse data
    m_crse_times[1] = crse_time[1]; // time of "new" coarse data

    m_dt_crse = crse_time[1] - crse_time[0];
}
