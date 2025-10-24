#ifndef NOAHMP_H
#define NOAHMP_H

#include <string>
#include <vector>
#include <memory>

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_Vector.H>
#include <AMReX_MFIter.H>

#include <ERF_NullSurf.H>
#include <ERF_Constants.H>
#include <ERF_IndexDefines.H>
#include <ERF_DataStruct.H>

// External include from the noahmp library
#include <NoahmpIO.H>

namespace LsmData_NOAHMP {
   enum {
         // 2D RRTMGP coupling vars
         t_sfc = 0,
         sfc_emis,
         sfc_alb_dir_vis,
         sfc_alb_dir_nir,
         sfc_alb_dif_vis,
         sfc_alb_dif_nir,
         cos_zenith_angle,
         sw_flux_dn,
         lw_flux_dn,
         NumVars
  };
}

namespace LsmFlux_NOAHMP {
   enum {
         // 2D SurfaceLayer flux vars
         t_flux = 0,
         q_flux,
         tau13,
         tau23,
         NumVars
  };
}

class NOAHMP : public NullSurf {

    using FabPtr = std::shared_ptr<amrex::MultiFab>;

public:
    // Constructor
    NOAHMP () {}

    // Destructor
    virtual ~NOAHMP () = default;

    // Set thermo and grid properties
    void
    Define (SolverChoice& /*sc*/) override
    {
        // NOTE: We should parse constants from sc here if needed,
    }

    // Initialize data structures
    void
    Init (const int& lev,
          const amrex::MultiFab& cons_in,
          const amrex::Geometry& geom,
          const amrex::Real& dt) override;

    // Wrapper to do all the updating
    void
    Advance_With_State (const int& lev,
                        amrex::MultiFab& cons_in,
                        amrex::MultiFab& xvel_in,
                        amrex::MultiFab& yvel_in,
                        amrex::MultiFab* hfx3_out,
                        amrex::MultiFab* qfx3_out,
                        const amrex::Real& dt,
                        const int& nstep) override;

    void
    Plot_Landfile (const int& nstep) override;

    // Get state vars from lsm class
    amrex::MultiFab*
    Lsm_Data_Ptr (const int& varIdx) override
    {
        int lsmIdx = LsmDataMap[varIdx];
        AMREX_ALWAYS_ASSERT(lsmIdx < NOAHMP::m_lsm_data_size && lsmIdx>=0);
        return lsm_fab_data[lsmIdx].get();
    }

    // Get flux vars from lsm class
    amrex::MultiFab*
    Lsm_Flux_Ptr (const int& varIdx) override
    {
        int lsmIdx = LsmFluxMap[varIdx];
        AMREX_ALWAYS_ASSERT(lsmIdx < NOAHMP::m_lsm_flux_size && lsmIdx>=0);
        return lsm_fab_flux[lsmIdx].get();
    }

    // Get lsm geometry
    amrex::Geometry
    Lsm_Geom ( ) override { return m_lsm_geom; }

    // Get number of vars lsm class contains
    int
    Lsm_Data_Size () override { return NOAHMP::m_lsm_data_size; }

    // Get number of fluxes lsm class contains
    int
    Lsm_Flux_Size () override { return NOAHMP::m_lsm_flux_size; }

    // Get variable names
    std::string
    Lsm_DataName (const int& varIdx) override
    {
        int lsmIdx = LsmDataMap[varIdx];
        AMREX_ALWAYS_ASSERT(lsmIdx < NOAHMP::m_lsm_data_size && lsmIdx>=0);
        return LsmDataName[lsmIdx];
    }

    // Get variable index from name
    int
    Lsm_DataIndex (std::string varname) override
    {
        int varIdx = -1;
        std::string lc_varname = amrex::toLower(varname);
        for (int idx(0); idx<LsmData_NOAHMP::NumVars; ++idx) {
            if (lc_varname == amrex::toLower(LsmDataName[idx])) {
                varIdx = idx;
            }
        }
        return varIdx;
    }

    // Get flux variable names
    std::string
    Lsm_FluxName (const int& varIdx) override
    {
        int lsmIdx = LsmFluxMap[varIdx];
        AMREX_ALWAYS_ASSERT(lsmIdx < NOAHMP::m_lsm_flux_size && lsmIdx>=0);
        return LsmFluxName[lsmIdx];
    }

    // Get flux variable index from name
    int
    Lsm_FluxIndex (std::string varname) override
    {
        int varIdx = -1;
        std::string lc_varname = amrex::toLower(varname);
        for (int idx(0); idx<LsmFlux_NOAHMP::NumVars; ++idx) {
            if (lc_varname == amrex::toLower(LsmFluxName[idx])) {
                varIdx = idx;
            }
        }
        return varIdx;
    }

private:
    // number of lsm data variables
    int m_lsm_data_size = LsmData_NOAHMP::NumVars;

    // number of lsm flux variables
    int m_lsm_flux_size = LsmFlux_NOAHMP::NumVars;

    // LsmData map (state indices -> LsmVar enum)
    amrex::Vector<int> LsmDataMap;

    // LsmFlux map (state indices -> LsmFlux enum)
    amrex::Vector<int> LsmFluxMap;

    // Lsm data names
    amrex::Vector<std::string> LsmDataName;

    // Lsm flux names
    amrex::Vector<std::string> LsmFluxName;

    // geometry for atmosphere
    amrex::Geometry m_geom;

    // geometry for lsm
    amrex::Geometry m_lsm_geom;

    // timestep
    amrex::Real m_dt;

    // domain klo-1 or lsm khi
    int khi_lsm;

    // Number of grid points in z
    int m_nz_lsm = 1;

    // Size of grid spacing in z
    amrex::Real m_dz_lsm = 1.0;

    // LSM data variables
    amrex::Array<FabPtr, LsmData_NOAHMP::NumVars> lsm_fab_data;

    // LSM data variables
    amrex::Array<FabPtr, LsmFlux_NOAHMP::NumVars> lsm_fab_flux;

    // Vector to store NoahmpIO objects at a level
    NoahmpIO_vector noahmpio_vect;

    // Plot frequency
    int m_plot_int_1 = -1;
};
#endif
