#ifndef MM5_H
#define MM5_H

#include <string>
#include <vector>
#include <memory>

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFabUtil.H>

#include <NullSurf.H>
#include <ERF_Constants.H>
#include <IndexDefines.H>
#include <DataStruct.H>

namespace LsmVar_MM5 {
   enum {
      // independent variables
      theta = 0,
      NumVars
  };
}

/* See Chen & Dudhia (2001) https://doi.org/10.1175/1520-0493(2001)129<0569:CAALSH>2.0.CO;2 */
class MM5 : public NullSurf {

    using FabPtr = std::shared_ptr<amrex::MultiFab>;

public:
    // Constructor
    MM5 () {}

    // Destructor
    virtual ~MM5 () = default;

    // Set thermo and grid properties
    void
    Define (SolverChoice& /*sc*/) override
    {
        // NOTE: We should parse things from sc here,
        //       but they are hard coded because this
        //       is a demonstration for now.
    }

    // Initialize data structures
    void
    Init (const amrex::MultiFab& cons_in,
          const amrex::Geometry& geom,
          const amrex::Real& dt) override;

    // Wrapper to do all the updating
    void
    Advance (const amrex::Real& dt) override
    {
        m_dt = dt;
        this->ComputeFluxes();
        this->AdvanceMM5();
        this->ComputeTsurf();
    }

    // Compute surface temperature
    void
    ComputeTsurf ();

    // Compute diffusive fluxes
    void
    ComputeFluxes ();

    // Advance the lsm state vars
    void
    AdvanceMM5 ();

    // Get state vars from lsm class
    amrex::MultiFab*
    Lsm_Data_Ptr (const int& varIdx) override
    {
        int lsmIdx = LsmVarMap[varIdx];
        AMREX_ALWAYS_ASSERT(lsmIdx < MM5::m_lsm_size && lsmIdx>=0);
        return lsm_fab_vars[lsmIdx].get();
    }

    // Get flux vars from lsm class
    amrex::MultiFab*
    Lsm_Flux_Ptr (const int& varIdx) override
    {
        int lsmIdx = LsmVarMap[varIdx];
        AMREX_ALWAYS_ASSERT(lsmIdx < MM5::m_lsm_size && lsmIdx>=0);
        return lsm_fab_flux[lsmIdx].get();
    }

    // Get lsm geometry
    amrex::Geometry
    Lsm_Geom ( ) override { return m_lsm_geom; }

    // Get number of vars lsm class contains
    int
    Lsm_Data_Size () override { return MM5::m_lsm_size; }

    // Get variable names
    std::string
    Lsm_VarName (const int& varIdx) override
    {
        int lsmIdx = LsmVarMap[varIdx];
        AMREX_ALWAYS_ASSERT(lsmIdx < MM5::m_lsm_size && lsmIdx>=0);
        return LsmVarName[lsmIdx];
    }

private:
    // number of lsm variables (theta)
    int m_lsm_size = 1;

    // LsmVar map (state indices -> LsmVar enum)
    amrex::Vector<int> LsmVarMap;

    // Lsm varnames
    amrex::Vector<std::string> LsmVarName;

    // geometry for atmosphere
    amrex::Geometry m_geom;

    // geometry for lsm
    amrex::Geometry m_lsm_geom;

    // timestep
    amrex::Real m_dt;

    // domain klo-1 or lsm khi
    int khi_lsm;

    // constants
    // amrex::Real m_fac_cond;
    // amrex::Real m_fac_fus;
    // amrex::Real m_fac_sub;
    // amrex::Real m_gOcp;

    // independent variables
    amrex::Array<FabPtr, LsmVar_MM5::NumVars> lsm_fab_vars;

    // flux array for conjugate transfer
    amrex::Array<FabPtr, LsmVar_MM5::NumVars> lsm_fab_flux;

    // Vars that should be parsed
    // ==========================
    // Number of grid points in z
    int m_nz_lsm = 30;

    // Size of grid spacing in z
    amrex::Real m_dz_lsm = 0.1; // 3 [m] below surface

    // Specific heat
    amrex::Real m_cp_soil =  1.26e6; // [J/m^3 K]

    // Conductivity
    amrex::Real m_k_soil = 0.2; // dry soil [W/m K]

    // Theta Dirichlet value at lowest point below surface
    amrex::Real m_theta_dir = 300.0; // [K]

    // Thermal diffusivity
    amrex::Real m_d_soil = m_k_soil / m_cp_soil; // [m^2/s]
};
#endif
