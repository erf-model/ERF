#ifndef ERF_FILLPATCHER_H_
#define ERF_FILLPATCHER_H_

#include <AMReX_FillPatchUtil.H>

class ERFFillPatcher
{
public:

    ERFFillPatcher (amrex::BoxArray const& fba, amrex::DistributionMapping const& fdm,
                    amrex::Geometry const& fgeom,
                    amrex::BoxArray const& cba, amrex::DistributionMapping const& cdm,
                    amrex::Geometry const& cgeom,
                    int nghost, int ncomp, amrex::InterpBase* interp);

    void registerCoarseData (amrex::Vector<amrex::MultiFab const*> const& crse_data,
                             amrex::Vector<amrex::Real> const& crse_time);

    template <typename BC>
    void fillCoarseFineBoundary (amrex::MultiFab& mf, amrex::Real time,
                                 BC& cbc, amrex::Vector<amrex::BCRec> const& bcs);

private:

    amrex::BoxArray m_fba;
    amrex::BoxArray m_cba;
    amrex::DistributionMapping m_fdm;
    amrex::DistributionMapping m_cdm;
    amrex::Geometry m_fgeom;
    amrex::Geometry m_cgeom;
    int m_nghost;
    int m_ncomp;
    amrex::InterpBase* m_interp;
    amrex::IntVect m_ratio;
    amrex::Vector<amrex::MultiFab> m_cf_crse_data;
    amrex::MultiFab m_cf_fine_data;
};

template <typename BC>
void
ERFFillPatcher::fillCoarseFineBoundary (amrex::MultiFab& mf, amrex::Real time,
                                        BC& cbc, amrex::Vector<amrex::BCRec> const& bcs)
{
    cbc(m_cf_crse_data[0], 0, m_ncomp, amrex::IntVect(0), time, 0);
    amrex::FillPatchInterp(m_cf_fine_data, 0, m_cf_crse_data[0], 0, m_ncomp, amrex::IntVect(0),
                           m_cgeom, m_fgeom, m_fgeom.Domain(), m_ratio, m_interp, bcs, 0);
    mf.ParallelCopy(m_cf_fine_data, 0, 0, m_ncomp, amrex::IntVect{0}, mf.nGrowVect());
}

#endif
