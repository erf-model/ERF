#ifndef ERF_RADIATION_SIMPLE_H
#define ERF_RADIATION_SIMPLE_H

#include <memory>

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFabUtil.H>

#include <ERF_RadiationInterface.H>
#include <ERF_DataStruct.H>
#include <ERF_EOS.H>

class RadiationSimple : public IRadiation
{
public:
    RadiationSimple (const int& /*lev*/,
                     const SolverChoice& sc) {
        // Every moist species we read is taken from the configured moisture-component
        // mapping, so each scheme is indexed by its own layout rather than by an
        // enumerated model list: every scheme carrying cloud ice (SAM, Morrison, WSM6,
        // WDM6, ...) is picked up, and the warm-rain schemes are not.
        m_qv_comp = sc.moisture_indices.qv;
        m_qc_comp = sc.moisture_indices.qc;
        m_qi_comp = sc.moisture_indices.qi;

        m_have_qv = (m_qv_comp >= 0);
        m_have_qc = (m_qc_comp >= 0);
        m_ice     = (m_qi_comp >= 0);
    }

    ~RadiationSimple () = default;

    virtual void Init (const amrex::Geometry& geom,
                       const amrex::BoxArray& ba,
                       amrex::MultiFab* cons_in) override;

    virtual void Run (int& level,
                      int& step,
                      double& time,
                      const double& dt,
                      const amrex::BoxArray& ba,
                      amrex::Geometry& geom,
                      amrex::MultiFab* cons_in,
                      amrex::iMultiFab* lmask,
                      amrex::MultiFab* t_surf,
                      const amrex::Vector<const amrex::MultiFab*>& lsm_input_ptrs,
                      const amrex::Vector<amrex::MultiFab*>& lsm_output_ptrs,
                      amrex::MultiFab* qheating_rates,
                      amrex::MultiFab* rad_fluxes,
                      amrex::MultiFab* z_phys,
                      amrex::MultiFab* z_sfc,
                      amrex::MultiFab* lat,
                      amrex::MultiFab* lon,
                      const bool force_update,
                      amrex::Real& declin,
                      amrex::Real& calday) override;

    virtual void WriteDataLog (const double &time) override;

private:
    bool m_have_qv = false;
    bool m_have_qc = false;
    bool m_ice     = false;
    int  m_qv_comp = -1;
    int  m_qc_comp = -1;
    int  m_qi_comp = -1;

    amrex::Geometry m_geom;
    amrex::BoxArray m_ba;
    bool m_partial_columns = false;

    amrex::MultiFab deltaq, flux;

    // Holds output for the datalog: cell-centered profiles, and the net flux on the surface
    // and top-of-atmosphere faces (one cell thick in z)
    amrex::MultiFab datalog_mf;
    amrex::MultiFab datalog_bdy_mf;

    // Layout of datalog_mf: height above ground, net lw flux (averaged from the faces) and
    // lw heating rate, all at cell centers
    static constexpr int dl_z_comp       = 0;
    static constexpr int dl_flux_comp    = 1;
    static constexpr int dl_heating_comp = 2;
    static constexpr int n_datalog_comp  = 3;

    // Layout of datalog_bdy_mf: net lw flux on the surface and top faces, then the top-face height
    static constexpr int bdy_sfc_comp       = 0;
    static constexpr int bdy_toa_comp       = 1;
    static constexpr int bdy_ztop_comp      = 2;
    static constexpr int n_datalog_bdy_comp = 3;
};

#endif
