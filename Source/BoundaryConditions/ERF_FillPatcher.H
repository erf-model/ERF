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
                    int nghost, int ncomp, amrex::InterpBase* interp);

    void registerCoarseData (amrex::Vector<amrex::MultiFab const*> const& crse_data,
                             amrex::Vector<amrex::Real> const& crse_time);

    template <typename BC>
    void fill(amrex::MultiFab& mf, amrex::Real time,
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
                      BC& cbc, amrex::Vector<amrex::BCRec> const& bcs)
{
    cbc(m_cf_crse_data[0], 0, m_ncomp, amrex::IntVect(0), time, 0);

    amrex::MultiFab crse_data_time_interp(m_cf_crse_data[0].boxArray(), m_cf_crse_data[0].DistributionMap(),
                                          m_ncomp, amrex::IntVect{0});

    amrex::Real fac_new = (time - m_crse_times[0]) / m_dt_crse;
    amrex::Real fac_old = 1.0 - fac_new;

    // Time interpolate the coarse data
    amrex::MultiFab::LinComb(crse_data_time_interp,fac_old, m_cf_crse_data[0], 0,
                                                   fac_new, m_cf_crse_data[1], 0,
                             0, m_ncomp, 0);


    // Spatially interpolate the time-interpolated coarse data
    amrex::FillPatchInterp(m_cf_fine_data, 0, crse_data_time_interp, 0, m_ncomp, amrex::IntVect(0),
                           m_cgeom, m_fgeom, m_fgeom.Domain(), m_ratio, m_interp, bcs, 0);

    // Copy the spatially and temporally interpolated data back into the mf
    mf.ParallelCopy(m_cf_fine_data, 0, 0, m_ncomp, amrex::IntVect{0}, mf.nGrowVect());
}
#endif
