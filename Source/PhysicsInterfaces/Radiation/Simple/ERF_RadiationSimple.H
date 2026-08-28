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
                      amrex::Vector<amrex::MultiFab*>& lsm_input_ptrs,
                      amrex::Vector<amrex::MultiFab*>& lsm_output_ptrs,
                      amrex::MultiFab* qheating_rates,
                      amrex::MultiFab* rad_fluxes,
                      amrex::MultiFab* z_phys,
                      amrex::MultiFab* lat,
                      amrex::MultiFab* lon,
                      const bool updated_lsm) override;

    virtual void WriteDataLog (const double &time) override;

    const amrex::MultiFab *get_radlwdn () { return radlwdn.get(); }
    const amrex::MultiFab *get_radqrlw () { return radqrlw.get(); }

private:
    bool m_have_qv = false;
    bool m_have_qc = false;
    bool m_ice     = false;
    int  m_qv_comp = -1;
    int  m_qc_comp = -1;
    int  m_qi_comp = -1;

    amrex::Geometry m_geom;
    amrex::BoxArray m_ba;

    amrex::MultiFab deltaq, flux;
    std::unique_ptr<amrex::MultiFab> radlwdn;
    std::unique_ptr<amrex::MultiFab> radqrlw;

};

#endif
