#include <ERF_FillPatcher.H>

using namespace amrex;

ERFFillPatcher::ERFFillPatcher (BoxArray const& fba, DistributionMapping const& fdm,
                                Geometry const& fgeom,
                                BoxArray const& cba, DistributionMapping const& cdm,
                                Geometry const& cgeom,
                                int nghost, int ncomp, InterpBase* interp)
    : m_fba(fba),
      m_cba(cba),
      m_fdm(fdm),
      m_cdm(cdm),
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
    m_cf_fine_data.define(cf_fba, cf_dm, m_ncomp, 0);
    m_cf_crse_data.emplace_back(BoxArray(std::move(cbl)), cf_dm, m_ncomp, 0);
}

void ERFFillPatcher::registerCoarseData (Vector<MultiFab const*> const& crse_data,
                                         Vector<amrex::Real> const& /*crse_time*/)
{
    AMREX_ALWAYS_ASSERT(crse_data.size() == 1);
    m_cf_crse_data[0].ParallelCopy(*crse_data[0], m_cgeom.periodicity());
}
