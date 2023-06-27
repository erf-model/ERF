#include <ERF_FillPatcher.H>
#include <utility>

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
                                int nghost, int ncomp, InterpBase* interp)
    : m_fba(fba),
      m_cba(std::move(cba)),
      m_fdm(std::move(fdm)),
      m_cdm(std::move(cdm)),
      m_fgeom(fgeom),
      m_cgeom(cgeom),
      m_nghost(nghost),
      m_ncomp(ncomp),
      m_interp(interp)
{
    AMREX_ALWAYS_ASSERT(nghost < 0);
    Box const fine_valid_box = amrex::grow(fba.minimalBox(), IntVect(nghost,nghost,0));

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        m_ratio[idim] = m_fgeom.Domain().length(idim) / m_cgeom.Domain().length(idim);
    }

    BoxList bl;
    const int nfboxes = fba.size();
    for (int ibox = 0; ibox < nfboxes; ++ibox) {
        Box vbx = fba[ibox];
        BoxList const& bndry = amrex::boxDiff(vbx, fine_valid_box);
        if (bndry.isNotEmpty()) {
            bl.join(bndry);
        }
    }

    BoxList cbl;
    cbl.reserve(bl.size());
    for (auto const& b : bl) {
        cbl.push_back(interp->CoarseBox(b, m_ratio));
    }

    BoxArray cf_fba(std::move(bl));
    DistributionMapping cf_dm(cf_fba);

    // This will be used as a temporary to hold the time-interpolated state
    m_cf_fine_data.define(cf_fba, cf_dm, m_ncomp, 0);

    // These will hold the coarse data on the m_cf grids
    m_cf_crse_data.emplace_back(BoxArray(std::move(cbl)), cf_dm, m_ncomp, 0);
    m_cf_crse_data.emplace_back(BoxArray(std::move(cbl)), cf_dm, m_ncomp, 0);
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
    AMREX_ALWAYS_ASSERT(crse_time[0] != crse_time[1]);

    m_cf_crse_data[0].ParallelCopy(*crse_data[0], m_cgeom.periodicity()); // old data
    m_cf_crse_data[1].ParallelCopy(*crse_data[1], m_cgeom.periodicity()); // new data

    m_crse_times.resize(2);
    m_crse_times[0] = crse_time[0]; // time of "old" coarse data
    m_crse_times[1] = crse_time[1]; // time of "new" coarse data

    m_dt_crse = crse_time[1] - crse_time[0];
}
