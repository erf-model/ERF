#ifndef ERF_FILLPATCHER_H_
#define ERF_FILLPATCHER_H_

#include <AMReX_FillPatchUtil.H>

class ERFFillPatcher
{
public:

    ERFFillPatcher (amrex::BoxArray const& fba, amrex::DistributionMapping  fdm,
                    amrex::Geometry const& fgeom,
                    amrex::BoxArray  cba, amrex::DistributionMapping  cdm,
                    amrex::Geometry const& cgeom,
                    int nghost, int nghost_subset, int ncomp, amrex::InterpBase* interp);

    ~ERFFillPatcher ( )
    {
        delete m_cf_fine_data;
        delete m_cf_fine_subset_data;
        delete m_cf_crse_data[0];
        delete m_cf_crse_data[1];
    }

    void Define (amrex::BoxArray const& fba, amrex::DistributionMapping  fdm,
                 amrex::Geometry const& fgeom,
                 amrex::BoxArray  cba, amrex::DistributionMapping  cdm,
                 amrex::Geometry const& cgeom,
                 int nghost, int nghost_subset, int ncomp, amrex::InterpBase* interp);

    void registerCoarseData (amrex::Vector<amrex::MultiFab const*> const& crse_data,
                             amrex::Vector<amrex::Real> const& crse_time);

    template <typename BC>
    void fill (amrex::MultiFab& mf, amrex::Real time,
              BC& cbc, amrex::Vector<amrex::BCRec> const& bcs, bool fill_subset=false);

private:

    amrex::BoxArray m_fba;
    amrex::BoxArray m_cba;
    amrex::DistributionMapping m_fdm;
    amrex::DistributionMapping m_cdm;
    amrex::Geometry m_fgeom;
    amrex::Geometry m_cgeom;
    int m_nghost;
    int m_nghost_subset;
    int m_ncomp;
    amrex::InterpBase* m_interp;
    amrex::IntVect m_ratio;
    amrex::Vector<amrex::MultiFab*> m_cf_crse_data;
    amrex::MultiFab* m_cf_fine_data;
    amrex::MultiFab* m_cf_fine_subset_data;
    amrex::Vector<amrex::Real> m_crse_times;
    amrex::Real m_dt_crse;
};

/*
 * Fill data at the coarse-fine boundary
 *
 * @param[out] mf    MultiFab to be filled
 * @param[in]  time  Time at which to fill data
 * @param[in]  cbc   Coarse boundary condition
 * @param[in]  bcs   Vector of boundary conditions
 */
template <typename BC>
void
ERFFillPatcher::fill (amrex::MultiFab& mf, amrex::Real time,
                      BC& cbc, amrex::Vector<amrex::BCRec> const& bcs, bool fill_subset)
{
    constexpr amrex::Real eps = std::numeric_limits<float>::epsilon();
    AMREX_ALWAYS_ASSERT((time >= m_crse_times[0]-eps) && (time <= m_crse_times[1]+eps));

    // Time interpolation factors
    amrex::Real fac_new = (time - m_crse_times[0]) / m_dt_crse;
    amrex::Real fac_old = 1.0 - fac_new;

    // Boundary condition operator
    cbc(*(m_cf_crse_data[0]), 0, m_ncomp, amrex::IntVect(0), time, 0);

    // Coarse MF to hold time interpolated data
    amrex::MultiFab crse_data_time_interp(m_cf_crse_data[0]->boxArray(), m_cf_crse_data[0]->DistributionMap(),
                                          m_ncomp, amrex::IntVect{0});

    // Time interpolate the coarse data
    amrex::MultiFab::LinComb(crse_data_time_interp,
                             fac_old, *(m_cf_crse_data[0]), 0,
                             fac_new, *(m_cf_crse_data[1]), 0,
                             0, m_ncomp, amrex::IntVect{0});

    // Ensure fine domain box is correct index type
    amrex::Box fdest_dom = amrex::convert(m_fgeom.Domain(),m_cf_fine_data->boxArray().ixType());

    // Spatially interpolate the time-interpolated coarse data
    amrex::FillPatchInterp(*m_cf_fine_data, 0, crse_data_time_interp, 0, m_ncomp, amrex::IntVect(0),
                           m_cgeom, m_fgeom, fdest_dom, m_ratio, m_interp, bcs, 0);

    // Fill whole region or subset?
    if (fill_subset) {
        //amrex::MultiFab::Copy(m_cf_fine_subset_data, m_cf_fine_data, 0, 0, m_ncomp, amrex::IntVect{0});
        m_cf_fine_subset_data->ParallelCopy(*m_cf_fine_data, 0, 0, m_ncomp, amrex::IntVect{0}, amrex::IntVect{0});
        mf.ParallelCopy(*m_cf_fine_subset_data, 0, 0, m_ncomp, amrex::IntVect{0}, amrex::IntVect{0});
    } else {
        mf.ParallelCopy(*m_cf_fine_data, 0, 0, m_ncomp, amrex::IntVect{0}, amrex::IntVect{0});
    }
}
#endif
