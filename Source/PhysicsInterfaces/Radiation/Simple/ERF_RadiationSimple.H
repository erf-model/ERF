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
    RadiationSimple(const int& /*lev*/,
                    const SolverChoice& sc) {
        // Check if we have a valid moisture model
        if (sc.moisture_type != MoistureType::None) { m_moist = true; }

        // Check if we have a moisture model with ice
        if (sc.moisture_type == MoistureType::SAM ||
            sc.moisture_type == MoistureType::Kessler ||
            sc.moisture_type == MoistureType::Morrison) { m_ice = true; }
    }

    ~RadiationSimple() = default;

    virtual void Init(const amrex::Geometry& geom,
                      const amrex::BoxArray& ba,
                      amrex::MultiFab* cons_in) override;

    virtual void Run(int& level,
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

    virtual void WriteDataLog(const double &time) override;

    const amrex::MultiFab *get_radlwdn() { return radlwdn.get(); }
    const amrex::MultiFab *get_radqrlw() { return radqrlw.get(); }

private:
    bool m_moist = false;
    bool m_ice = false;

    amrex::Geometry m_geom;
    amrex::BoxArray m_ba;

    amrex::MultiFab deltaq, flux;
    std::unique_ptr<amrex::MultiFab> radlwdn;
    std::unique_ptr<amrex::MultiFab> radqrlw;

};

#endif
