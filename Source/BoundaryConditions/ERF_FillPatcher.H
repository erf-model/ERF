#ifndef ERF_FILLPATCHER_H_
#define ERF_FILLPATCHER_H_

#include <AMReX_FillPatchUtil.H>
#include <AMReX_Interp_C.H>
#include <AMReX_MFInterp_C.H>

class ERFFillPatcher
{
public:

    /**
     * Construct a coarse-fine fill patcher.
     *
     * @param[in] fba         fine-level BoxArray
     * @param[in] fdm         fine-level DistributionMapping
     * @param[in] fgeom       fine-level Geometry
     * @param[in] cba         coarse-level BoxArray
     * @param[in] cdm         coarse-level DistributionMapping
     * @param[in] cgeom       coarse-level Geometry
     * @param[in] nghost      number of ghost cells in the fill region
     * @param[in] nghost_set  number of ghost cells in the set region
     * @param[in] ncomp       number of components
     * @param[in] interp      interpolation operator
     */
    ERFFillPatcher (amrex::BoxArray const& fba, amrex::DistributionMapping const& fdm,
                    amrex::Geometry const& fgeom,
                    amrex::BoxArray const& cba, amrex::DistributionMapping const& cdm,
                    amrex::Geometry const& cgeom,
                    int nghost, int nghost_set, int ncomp, amrex::InterpBase* interp);

    /**
     * Define the coarse-fine fill patcher.
     *
     * @param[in] fba         fine-level BoxArray
     * @param[in] fdm         fine-level DistributionMapping
     * @param[in] fgeom       fine-level Geometry
     * @param[in] cba         coarse-level BoxArray
     * @param[in] cdm         coarse-level DistributionMapping
     * @param[in] cgeom       coarse-level Geometry
     * @param[in] nghost      number of ghost cells in the fill region
     * @param[in] nghost_set  number of ghost cells in the set region
     * @param[in] ncomp       number of components
     * @param[in] interp      interpolation operator
     */
    void Define (amrex::BoxArray const& fba, amrex::DistributionMapping const& fdm,
                 amrex::Geometry const& fgeom,
                 amrex::BoxArray const& cba, amrex::DistributionMapping const& cdm,
                 amrex::Geometry const& cgeom,
                 int nghost, int nghost_set, int ncomp,
                 amrex::InterpBase* interp);

    /**
     * Build the coarse-fine mask.
     *
     * @param[in] fba        fine-level BoxArray
     * @param[in] nghost     number of ghost cells in the fill region
     * @param[in] nghost_set number of ghost cells in the set region
     */
    void BuildMask (amrex::BoxArray const& fba, int nghost, int nghost_set);

    /**
     * Register coarse data bracketing the fill time.
     *
     * @param[in] crse_data coarse MultiFabs
     * @param[in] crse_time times associated with crse_data
     */
    void RegisterCoarseData (amrex::Vector<amrex::MultiFab const*> const& crse_data,
                             amrex::Vector<double> const& crse_time);

    /**
     * Interpolate coarse data to fine face data.
     *
     * @param[in,out] fine     fine-level data to fill
     * @param[in]     crse     coarse-level source data
     * @param[in]     mask_val mask value selecting cells to fill
     */
    void InterpFace (amrex::MultiFab& fine,
                     amrex::MultiFab const& crse,
                     int mask_val);

    /**
     * Interpolate coarse data to fine cell-centered data.
     *
     * @param[in,out] fine     fine-level data to fill
     * @param[in]     crse     coarse-level source data
     * @param[in]     bcr      boundary-condition records
     * @param[in]     mask_val mask value selecting cells to fill
     */
    void InterpCell (amrex::MultiFab& fine,
                     amrex::MultiFab const& crse,
                     amrex::Vector<amrex::BCRec> const& bcr,
                     int mask_val);

    /**
     * Return the mask value for the set region.
     */
    int GetSetMaskVal () { return m_set_mask; }

    /**
     * Return the mask value for the relax region.
     */
    int GetRelaxMaskVal () { return m_relax_mask; }

    /**
     * Return the coarse-fine mask.
     */
    amrex::iMultiFab* GetMask () { return m_cf_mask.get(); }

    /**
     * Fill fine data in the set region.
     *
     * @param[in,out] mf   MultiFab to be filled
     * @param[in]     time time at which to fill data
     * @param[in]     cbc  coarse boundary condition
     * @param[in]     bcs  vector of boundary conditions
     */
    template <typename BC>
    void FillSet (amrex::MultiFab& mf, double time,
                  BC& cbc, amrex::Vector<amrex::BCRec> const& bcs);

    /**
     * Fill fine data in the relax region.
     *
     * @param[in,out] mf   MultiFab to be filled
     * @param[in]     time time at which to fill data
     * @param[in]     cbc  coarse boundary condition
     * @param[in]     bcs  vector of boundary conditions
     */
    template <typename BC>
    void FillRelax (amrex::MultiFab& mf, double time,
                    BC& cbc, amrex::Vector<amrex::BCRec> const& bcs);

    /**
     * Fill fine data in the selected mask region.
     *
     * @param[in,out] mf       MultiFab to be filled
     * @param[in]     time     time at which to fill data
     * @param[in]     cbc      coarse boundary condition
     * @param[in]     bcs      vector of boundary conditions
     * @param[in]     mask_val mask value selecting cells to fill
     */
    template <typename BC>
    void Fill (amrex::MultiFab& mf, double time,
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
    amrex::Vector<double> m_crse_times;
    double m_dt_crse;
    int m_set_mask{2};
    int m_relax_mask{1};
};

/**
 * Fill fine data in the set region
 *
 * @param[in,out] mf   MultiFab to be filled
 * @param[in]     time Time at which to fill data
 * @param[in]     cbc  Coarse boundary condition
 * @param[in]     bcs  Vector of boundary conditions
 */
template <typename BC>
void
ERFFillPatcher::FillSet (amrex::MultiFab& mf, double time,
                         BC& cbc, amrex::Vector<amrex::BCRec> const& bcs)
{
    Fill(mf,time,cbc,bcs,m_set_mask);
}

/**
 * Fill fine data in the relax region
 *
 * @param[in,out] mf   MultiFab to be filled
 * @param[in]     time Time at which to fill data
 * @param[in]     cbc  Coarse boundary condition
 * @param[in]     bcs  Vector of boundary conditions
 */
template <typename BC>
void
ERFFillPatcher::FillRelax (amrex::MultiFab& mf, double time,
                           BC& cbc, amrex::Vector<amrex::BCRec> const& bcs)
{
    Fill(mf,time,cbc,bcs,m_relax_mask);
}

/**
 * Fill fine data in the relax region
 *
 * @param[in,out] mf       MultiFab to be filled
 * @param[in]     time_d   Time at which to fill data
 * @param[in]     cbc      Coarse boundary condition
 * @param[in]     bcs      Vector of boundary conditions
 * @param[in]     mask_val Mask value selecting cells to fill
 */
template <typename BC>
void
ERFFillPatcher::Fill (amrex::MultiFab& mf, double time_d,
                      BC& cbc, amrex::Vector<amrex::BCRec> const& bcs, int mask_val)
{
    constexpr double eps = std::numeric_limits<float>::epsilon();

    AMREX_ALWAYS_ASSERT((time_d >= m_crse_times[0]-eps) && (time_d <= m_crse_times[1]+eps));

    amrex::Real time = static_cast<amrex::Real>(time_d);

    // Time interpolation factors
    amrex::Real fac_new = amrex::Real(0);
    amrex::Real fac_old = amrex::Real(1);
    if (amrex::Math::abs(m_dt_crse) > eps) {
        fac_new = static_cast<amrex::Real>((time_d - m_crse_times[0]) / m_dt_crse);
        fac_old = amrex::Real(1) - fac_new;
    }

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
