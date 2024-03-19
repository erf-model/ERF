#ifndef ERF_FILLPATCHER_H_
#define ERF_FILLPATCHER_H_

#include <AMReX_FillPatchUtil.H>
#include <AMReX_Interp_C.H>
#include <AMReX_MFInterp_C.H>

class ERFFillPatcher
{
public:

    ERFFillPatcher (amrex::BoxArray const& fba, amrex::DistributionMapping const& fdm,
                    amrex::Geometry const& fgeom,
                    amrex::BoxArray const& cba, amrex::DistributionMapping const& cdm,
                    amrex::Geometry const& cgeom,
                    int nghost, int nghost_set, int ncomp, amrex::InterpBase* interp);

    void Define (amrex::BoxArray const& fba, amrex::DistributionMapping const& fdm,
                 amrex::Geometry const& fgeom,
                 amrex::BoxArray const& cba, amrex::DistributionMapping const& cdm,
                 amrex::Geometry const& cgeom,
                 int nghost, int nghost_set, int ncomp,
                 amrex::InterpBase* interp);

    void BuildMask (amrex::BoxArray const& fba, int nghost, int nghost_set);

    void RegisterCoarseData (amrex::Vector<amrex::MultiFab const*> const& crse_data,
                             amrex::Vector<amrex::Real> const& crse_time);

    void InterpFace (amrex::MultiFab& fine,
                     amrex::MultiFab const& crse,
                     int mask_val);

    void InterpCell (amrex::MultiFab& fine,
                     amrex::MultiFab const& crse,
                     amrex::Vector<amrex::BCRec> const& bcr,
                     int mask_val);

    int GetSetMaskVal () { return m_set_mask; }

    int GetRelaxMaskVal () { return m_relax_mask; }

    amrex::iMultiFab* GetMask () { return m_cf_mask.get(); }

    template <typename BC>
    void FillSet (amrex::MultiFab& mf, amrex::Real time,
                  BC& cbc, amrex::Vector<amrex::BCRec> const& bcs);

    template <typename BC>
    void FillRelax (amrex::MultiFab& mf, amrex::Real time,
                    BC& cbc, amrex::Vector<amrex::BCRec> const& bcs);

    template <typename BC>
    void Fill (amrex::MultiFab& mf, amrex::Real time,
               BC& cbc, amrex::Vector<amrex::BCRec> const& bcs, int mask_val);

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
    std::unique_ptr<amrex::MultiFab> m_cf_crse_data_old;
    std::unique_ptr<amrex::MultiFab> m_cf_crse_data_new;
    std::unique_ptr<amrex::iMultiFab> m_cf_mask;
    amrex::Vector<amrex::Real> m_crse_times;
    amrex::Real m_dt_crse;
    int m_set_mask{2};
    int m_relax_mask{1};
};

/*
 * Fill fine data in the set region
 *
 * @param[out] mf    MultiFab to be filled
 * @param[in]  time  Time at which to fill data
 * @param[in]  cbc   Coarse boundary condition
 * @param[in]  bcs   Vector of boundary conditions
 */
template <typename BC>
void
ERFFillPatcher::FillSet (amrex::MultiFab& mf, amrex::Real time,
                         BC& cbc, amrex::Vector<amrex::BCRec> const& bcs)
{
    Fill(mf,time,cbc,bcs,m_set_mask);
}

/*
 * Fill fine data in the relax region
 *
 * @param[out] mf    MultiFab to be filled
 * @param[in]  time  Time at which to fill data
 * @param[in]  cbc   Coarse boundary condition
 * @param[in]  bcs   Vector of boundary conditions
 */
template <typename BC>
void
ERFFillPatcher::FillRelax (amrex::MultiFab& mf, amrex::Real time,
                           BC& cbc, amrex::Vector<amrex::BCRec> const& bcs)
{
    Fill(mf,time,cbc,bcs,m_relax_mask);
}

/*
 * Fill fine data in the relax region
 *
 * @param[out] mf    MultiFab to be filled
 * @param[in]  time  Time at which to fill data
 * @param[in]  cbc   Coarse boundary condition
 * @param[in]  bcs   Vector of boundary conditions
 * @param[in]  mask_val Value to assign mask array
 */
template <typename BC>
void
ERFFillPatcher::Fill (amrex::MultiFab& mf, amrex::Real time,
                      BC& cbc, amrex::Vector<amrex::BCRec> const& bcs, int mask_val)
{
    constexpr amrex::Real eps = std::numeric_limits<float>::epsilon();

    AMREX_ALWAYS_ASSERT((time >= m_crse_times[0]-eps) && (time <= m_crse_times[1]+eps));

    // Time interpolation factors
    amrex::Real fac_new = (time - m_crse_times[0]) / m_dt_crse;
    amrex::Real fac_old = 1.0 - fac_new;

    // Boundary condition operator
    cbc(*(m_cf_crse_data_old), 0, m_ncomp, amrex::IntVect(0), time, 0);

    // Coarse MF to hold time interpolated data
    amrex::MultiFab crse_data_time_interp(m_cf_crse_data_old->boxArray(),
                                          m_cf_crse_data_old->DistributionMap(),
                                          m_ncomp, amrex::IntVect{0});

    // Time interpolate the coarse data
    amrex::MultiFab::LinComb(crse_data_time_interp,
                             fac_old, *(m_cf_crse_data_old), 0,
                             fac_new, *(m_cf_crse_data_new), 0,
                             0, m_ncomp, amrex::IntVect{0});

    // Call correct spatial interpolation type
    amrex::IndexType m_ixt = mf.boxArray().ixType();
    int ixt_sum = m_ixt[0]+m_ixt[1]+m_ixt[2];
    if (ixt_sum == 0) {
        InterpCell(mf,crse_data_time_interp,bcs,mask_val);
    } else if (ixt_sum == 1) {
        InterpFace(mf,crse_data_time_interp,mask_val);
    } else {
        amrex::Abort("ERF_FillPatcher only supports face linear and cell cons linear interp!");
    }
}
#endif
