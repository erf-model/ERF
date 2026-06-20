#ifndef NOAHMP_H
#define NOAHMP_H

#include <string>
#include <vector>
#include <memory>
#include <limits>

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
         sw_flux_dn_dir_vis,
         sw_flux_dn_dir_nir,
         sw_flux_dn_dif_vis,
         sw_flux_dn_dif_nir,
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


namespace NoahmpInputComp {
    enum {
        u_phy = 0,
        v_phy,
        t_phy,
        qv_curr,
        p8w,
        swdown,
        glw,
        coszen,
        NumComps
    };
}

namespace NoahmpOutputComp {
    enum {
        hfx = 0,
        lh,
        tau_ew,
        tau_ns,
        tsk,
        emiss,
        albsfcdir_vis,
        albsfcdir_nir,
        albsfcdif_vis,
        albsfcdif_nir,
        NumComps
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
                        const amrex::Real& elapsed_time,
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

    // Checkpoint/restore the NoahMP substep counter. Uses the class member
    // m_itimestep (valid on every rank) rather than noahmpio_vect[0] (absent on
    // land-free ranks).
    int
    Get_LSM_Step () const override
    {
        return m_itimestep;
    }

    void
    Set_LSM_Step (int step) override
    {
        m_itimestep = step;
        for (auto& noahmpio : noahmpio_vect) {
            noahmpio.itimestep = step;
        }
    }

    // Write/read the full NoahMP prognostic state (soil, snow incl. layers,
    // canopy, aquifer, albedo history, ...) to/from a NetCDF restart directory
    // so a restart reproduces a cold-start trajectory bitwise (issue #3255).
    // Each local block writes its tile into the global-domain file collectively.
    void
    Write_Lsm_Restart (const std::string& dir) const override
    {
        for (const NoahmpIO_type& noahmpio : noahmpio_vect) {
            const_cast<NoahmpIO_type&>(noahmpio).WriteRestart(dir);
        }
    }

    void
    Read_Lsm_Restart (const std::string& dir) override
    {
        for (NoahmpIO_type& noahmpio : noahmpio_vect) {
            noahmpio.ReadRestart(dir);
        }
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

    // Authoritative NoahMP substep counter, kept on the class so the firing
    // decision in Advance_With_State is identical on every rank (including
    // land-free ranks, whose noahmpio_vect is empty).
    int m_itimestep = 0;

    // NoahMP timestep (DTBL), broadcast to all ranks in Init so the firing
    // decision is well-defined on land-free ranks too.
    amrex::Real m_dtbl = std::numeric_limits<amrex::Real>::max();

    // domain klo-1 or lsm khi
    int khi_lsm;

    // Number of grid points in z
    int m_nz_lsm = 1;

    // Size of grid spacing in z
    amrex::Real m_dz_lsm = one;

    // LSM data variables
    amrex::Array<FabPtr, LsmData_NOAHMP::NumVars> lsm_fab_data;

    // LSM data variables
    amrex::Array<FabPtr, LsmFlux_NOAHMP::NumVars> lsm_fab_flux;

    // Vector to store NoahmpIO objects at a level
    NoahmpIO_vector noahmpio_vect;

    // Pinned memory buffers for host/device coupling to Noah-MP
    amrex::Vector<std::unique_ptr<amrex::FArrayBox>> noahmp_input_tmp;
    amrex::Vector<std::unique_ptr<amrex::FArrayBox>> noahmp_output_tmp;

    // Plot frequency
    int m_plot_int_1 = -1;
};
#endif
