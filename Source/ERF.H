#ifndef ERF_H_
#define ERF_H_

#include <string>
#include <limits>
#include <memory>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <AMReX_AmrCore.H>
#include <AMReX_BCRec.H>

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_YAFluxRegister.H>
#include <AMReX_ErrorList.H>

#ifdef AMREX_MEM_PROFILING
#include <AMReX_MemProfiler.H>
#endif

#include "prob_common.H"

#include <IndexDefines.H>
#include <DataStruct.H>
#include <InputSoundingData.H>
#include <InputSpongeData.H>
#include <ABLMost.H>
#include <Derive.H>
#include <ERF_ReadBndryPlanes.H>
#include <ERF_WriteBndryPlanes.H>
#include <ERF_MRI.H>
#include <ERF_PhysBCFunct.H>
#include <ERF_FillPatcher.H>

#ifdef ERF_USE_PARTICLES
#include "ParticleData.H"
#endif

#include "EulerianMicrophysics.H"
#include "LagrangianMicrophysics.H"
#include "LandSurface.H"

#ifdef ERF_USE_RRTMGP
#include "Radiation.H"
#endif

#include <iostream>

#ifdef AMREX_LAZY
#include <AMReX_Lazy.H>
#endif

#ifdef ERF_USE_MULTIBLOCK
class MultiBlockContainer;
#endif

/**
 * Enum of possible coarse/fine interpolation options
*/
namespace InterpType {
    enum {
        PCInterp = 0,
        NodeBilinear,
        CellConservativeLinear,
        CellBilinear,
        CellQuadratic,
        CellConservativeProtected,
        CellConservativeQuartic,
        FaceLinear,
        FaceConserativeLinear
    };
}

/**
 * Main class in ERF code, instantiated from main.cpp
*/

class ERF
    : public amrex::AmrCore
{
public:

    ////////////////
    // public member functions

    // constructor - reads in parameters from inputs file
    //             - sizes multilevel arrays and data structures
    ERF ();
    ~ERF () override;

    // Declare a default move constructor so we ensure the destructor is
    // not called when we return an object of this class by value
    ERF (ERF&&)  noexcept = delete;

    // Declare a default move assignment operator
    ERF& operator=(ERF&& other)  noexcept = delete;

    // Delete the copy constructor
    ERF (const ERF& other) = delete;
    //
    // Delete the copy assignment operator
    ERF& operator=(const ERF& other) = delete;

    // Advance solution to final time
    void Evolve ();

    // Tag cells for refinement
    void ErrorEst (int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow) override;

    // Initialize multilevel data
    void InitData ();

#ifdef ERF_USE_EB
    void WriteMyEBSurface ();
#endif

#ifdef ERF_USE_POISSON_SOLVE
    // Project the velocities to be divergence-free
    void project_velocities (                                        amrex::Vector<amrex::MultiFab >& vars);
    void project_velocities (int lev_min, int lev_max, amrex::Vector<amrex::Vector<amrex::MultiFab>>& vars);

    // Define the projection bc's based on the domain bc types
    amrex::Array<amrex::LinOpBCType,AMREX_SPACEDIM>
      get_projection_bc (amrex::Orientation::Side side) const noexcept;
    bool projection_has_dirichlet (amrex::Array<amrex::LinOpBCType,AMREX_SPACEDIM> bcs) const;

    // Project the velocities to be divergence-free with a thin body
    void project_velocities_tb (              amrex::Vector<amrex::MultiFab >& vars, const amrex::Real dt);
    void project_velocities_tb (amrex::Vector<amrex::Vector<amrex::MultiFab>>& vars, const amrex::Real dt);
#endif

    // Init (NOT restart or regrid)
    void init_only (int lev, amrex::Real time);

    // Restart
    void restart ();

    // Is it time to write a plotfile or checkpoint?
    bool writeNow (const amrex::Real cur_time, const amrex::Real dt, const int nstep, const int plot_int, const amrex::Real plot_per);

    // Called after every level 0 timestep
    void post_timestep (int nstep, amrex::Real time, amrex::Real dt_lev);

    // Diagnostics
    void sum_integrated_quantities (amrex::Real time);
    void write_1D_profiles (amrex::Real time);
    void write_1D_profiles_stag (amrex::Real time);

    // Fill the physical boundary conditions for cell-centered velocity (diagnostic only)
    void FillBdyCCVels (amrex::Vector<amrex::MultiFab>& mf_cc_vel);

    void sample_points (int lev, amrex::Real time, amrex::IntVect cell, amrex::MultiFab& mf);
    void sample_lines (int lev, amrex::Real time, amrex::IntVect cell, amrex::MultiFab& mf);

    void derive_diag_profiles (amrex::Gpu::HostVector<amrex::Real>& h_avg_u   , amrex::Gpu::HostVector<amrex::Real>& h_avg_v  ,
                               amrex::Gpu::HostVector<amrex::Real>& h_avg_w   , amrex::Gpu::HostVector<amrex::Real>& h_avg_rho,
                               amrex::Gpu::HostVector<amrex::Real>& h_avg_th  , amrex::Gpu::HostVector<amrex::Real>& h_avg_ksgs,
                               amrex::Gpu::HostVector<amrex::Real>& h_avg_qv  , amrex::Gpu::HostVector<amrex::Real>& h_avg_qc ,
                               amrex::Gpu::HostVector<amrex::Real>& h_avg_qr  ,
                               amrex::Gpu::HostVector<amrex::Real>& h_avg_wqv , amrex::Gpu::HostVector<amrex::Real>& h_avg_wqc ,
                               amrex::Gpu::HostVector<amrex::Real>& h_avg_wqr ,
                               amrex::Gpu::HostVector<amrex::Real>& h_avg_uu  , amrex::Gpu::HostVector<amrex::Real>& h_avg_uv ,
                               amrex::Gpu::HostVector<amrex::Real>& h_avg_uw,
                               amrex::Gpu::HostVector<amrex::Real>& h_avg_vv  , amrex::Gpu::HostVector<amrex::Real>& h_avg_vw ,
                               amrex::Gpu::HostVector<amrex::Real>& h_avg_ww,
                               amrex::Gpu::HostVector<amrex::Real>& h_avg_uth , amrex::Gpu::HostVector<amrex::Real>& h_avg_vth,
                               amrex::Gpu::HostVector<amrex::Real>& h_avg_wth, amrex::Gpu::HostVector<amrex::Real>& h_avg_thth,
                               amrex::Gpu::HostVector<amrex::Real>& h_avg_ku, amrex::Gpu::HostVector<amrex::Real>& h_avg_kv,
                               amrex::Gpu::HostVector<amrex::Real>& h_avg_kw,
                               amrex::Gpu::HostVector<amrex::Real>& h_avg_p,
                               amrex::Gpu::HostVector<amrex::Real>& h_avg_pu, amrex::Gpu::HostVector<amrex::Real>& h_avg_pv,
                               amrex::Gpu::HostVector<amrex::Real>& h_avg_pw);
    void derive_diag_profiles_stag (amrex::Gpu::HostVector<amrex::Real>& h_avg_u   , amrex::Gpu::HostVector<amrex::Real>& h_avg_v  ,
                                    amrex::Gpu::HostVector<amrex::Real>& h_avg_w   , amrex::Gpu::HostVector<amrex::Real>& h_avg_rho,
                                    amrex::Gpu::HostVector<amrex::Real>& h_avg_th  , amrex::Gpu::HostVector<amrex::Real>& h_avg_ksgs,
                                    amrex::Gpu::HostVector<amrex::Real>& h_avg_uu  , amrex::Gpu::HostVector<amrex::Real>& h_avg_uv ,
                                    amrex::Gpu::HostVector<amrex::Real>& h_avg_uw,
                                    amrex::Gpu::HostVector<amrex::Real>& h_avg_vv  , amrex::Gpu::HostVector<amrex::Real>& h_avg_vw ,
                                    amrex::Gpu::HostVector<amrex::Real>& h_avg_ww,
                                    amrex::Gpu::HostVector<amrex::Real>& h_avg_uth , amrex::Gpu::HostVector<amrex::Real>& h_avg_vth,
                                    amrex::Gpu::HostVector<amrex::Real>& h_avg_wth, amrex::Gpu::HostVector<amrex::Real>& h_avg_thth,
                                    amrex::Gpu::HostVector<amrex::Real>& h_avg_ku, amrex::Gpu::HostVector<amrex::Real>& h_avg_kv,
                                    amrex::Gpu::HostVector<amrex::Real>& h_avg_kw,
                                    amrex::Gpu::HostVector<amrex::Real>& h_avg_p,
                                    amrex::Gpu::HostVector<amrex::Real>& h_avg_pu, amrex::Gpu::HostVector<amrex::Real>& h_avg_pv,
                                    amrex::Gpu::HostVector<amrex::Real>& h_avg_pw);

    void derive_stress_profiles (amrex::Gpu::HostVector<amrex::Real>& h_avg_tau11, amrex::Gpu::HostVector<amrex::Real>& h_avg_tau12,
                                 amrex::Gpu::HostVector<amrex::Real>& h_avg_tau13, amrex::Gpu::HostVector<amrex::Real>& h_avg_tau22,
                                 amrex::Gpu::HostVector<amrex::Real>& h_avg_tau23, amrex::Gpu::HostVector<amrex::Real>& h_avg_tau33,
                                 amrex::Gpu::HostVector<amrex::Real>& h_avg_hfx3,  amrex::Gpu::HostVector<amrex::Real>& h_avg_diss);
    void derive_stress_profiles_stag (amrex::Gpu::HostVector<amrex::Real>& h_avg_tau11, amrex::Gpu::HostVector<amrex::Real>& h_avg_tau12,
                                      amrex::Gpu::HostVector<amrex::Real>& h_avg_tau13, amrex::Gpu::HostVector<amrex::Real>& h_avg_tau22,
                                      amrex::Gpu::HostVector<amrex::Real>& h_avg_tau23, amrex::Gpu::HostVector<amrex::Real>& h_avg_tau33,
                                      amrex::Gpu::HostVector<amrex::Real>& h_avg_hfx3,  amrex::Gpu::HostVector<amrex::Real>& h_avg_diss);

    // Perform the volume-weighted sum
    amrex::Real
    volWgtSumMF (int lev, const amrex::MultiFab& mf, int comp,
                 const amrex::MultiFab& mapfac,
                 bool local, bool finemask);

    // Decide if it is time to take an action
    static bool is_it_time_for_action (int nstep, amrex::Real time, amrex::Real dt,
                                       int action_interval, amrex::Real action_per);

    // Make a new level using provided BoxArray and DistributionMapping and
    // fill with interpolated coarse level data.
    // overrides the pure virtual function in AmrCore
    void MakeNewLevelFromCoarse (int lev, amrex::Real time, const amrex::BoxArray& ba,
                                 const amrex::DistributionMapping& dm) override;

    // Remake an existing level using provided BoxArray and DistributionMapping and
    // fill with existing fine and coarse data.
    // overrides the pure virtual function in AmrCore
    void RemakeLevel (int lev, amrex::Real time, const amrex::BoxArray& ba,
                      const amrex::DistributionMapping& dm) override;

    // Delete level data
    // overrides the pure virtual function in AmrCore
    void ClearLevel (int lev) override;

    // Make a new level from scratch using provided BoxArray and DistributionMapping.
    // Only used during initialization.
    // overrides the pure virtual function in AmrCore
    void MakeNewLevelFromScratch (int lev, amrex::Real time, const amrex::BoxArray& ba,
                                  const amrex::DistributionMapping& dm) override;

    // compute dt from CFL considerations
    amrex::Real estTimeStep (int lev, long& dt_fast_ratio) const;

    // Interface for advancing the data at one level by one "slow" timestep
    void advance_dycore (int level,
                         amrex::Vector<amrex::MultiFab>& state_old,
                         amrex::Vector<amrex::MultiFab>& state_new,
                         amrex::MultiFab& xvel_old,  amrex::MultiFab& yvel_old,  amrex::MultiFab& zvel_old,
                         amrex::MultiFab& xvel_new,  amrex::MultiFab& yvel_new,  amrex::MultiFab& zvel_new,
                         amrex::MultiFab& source,    amrex::MultiFab& xmom_src,
                         amrex::MultiFab& ymom_src,  amrex::MultiFab& zmom_src,
                         amrex::Geometry fine_geom,
                         amrex::Real dt, amrex::Real time);

    void advance_microphysics (int lev,
                               amrex::MultiFab& cons_in,
                               const amrex::Real& dt_advance,
                               const int& iteration,
                               const amrex::Real& time);

    void advance_lsm (int lev,
                      amrex::MultiFab& /*cons_in*/,
                      const amrex::Real& dt_advance);

#if defined(ERF_USE_RRTMGP)
    void advance_radiation (int lev,
                            amrex::MultiFab& cons_in,
                            const amrex::Real& dt_advance);
#endif

    amrex::MultiFab& build_fine_mask (int lev);

    void MakeHorizontalAverages ();
    void MakeDiagnosticAverage (amrex::Vector<amrex::Real>& h_havg, amrex::MultiFab& S, int n);
    void derive_upwp (amrex::Vector<amrex::Real>& h_havg);

    // write plotfile to disk
    void WritePlotFile  (int which, amrex::Vector<std::string> plot_var_names);

    void WriteMultiLevelPlotfileWithTerrain (const std::string &plotfilename,
                                             int nlevels,
                                             const amrex::Vector<const amrex::MultiFab*> &mf,
                                             const amrex::Vector<const amrex::MultiFab*> &mf_nd,
                                             const amrex::Vector<std::string> &varnames,
                                             amrex::Real time,
                                             const amrex::Vector<int> &level_steps,
                                             const std::string &versionName = "HyperCLaw-V1.1",
                                             const std::string &levelPrefix = "Level_",
                                             const std::string &mfPrefix = "Cell",
                                             const amrex::Vector<std::string>& extra_dirs = amrex::Vector<std::string>()) const;


    void WriteGenericPlotfileHeaderWithTerrain (std::ostream &HeaderFile,
                                                int nlevels,
                                                const amrex::Vector<amrex::BoxArray> &bArray,
                                                const amrex::Vector<std::string> &varnames,
                                                amrex::Real time,
                                                const amrex::Vector<int> &level_steps,
                                                const std::string &versionName,
                                                const std::string &levelPrefix,
                                                const std::string &mfPrefix) const;

    void erf_enforce_hse (int lev,
                          amrex::MultiFab& dens, amrex::MultiFab& pres, amrex::MultiFab& pi,
                          std::unique_ptr<amrex::MultiFab>& z_cc);

#ifdef ERF_USE_NETCDF
    //! Write a timestep to 1D vertical column output for coupling
    void writeToNCColumnFile (int lev,
                              const std::string& colfile_name, amrex::Real xloc, amrex::Real yloc,
                              amrex::Real time);
#endif //ERF_USE_NETCDF

    void init_from_input_sounding (int lev);

    void input_sponge (int lev);

    void init_from_hse (int lev);

#ifdef ERF_USE_MULTIBLOCK
    // constructor used when ERF is created by a multiblock driver
    // calls AmrCore constructor -> AmrMesh constructor
    ERF (const amrex::RealBox& rb, int max_level_in,
         const amrex::Vector<int>& n_cell_in, int coord,
         const amrex::Vector<amrex::IntVect>& ref_ratios,
         const amrex::Array<int,AMREX_SPACEDIM>& is_per,
         std::string prefix);

    // Advance a block specified number of time steps
    void Evolve_MB (int MBstep, int max_block_step);

    // get the current time values
    amrex::Real get_t_old() {return t_old[0];}
    amrex::Real get_t_new() {return t_new[0];}

    // Set parmparse prefix for MultiBlock
    void SetParmParsePrefix (std::string name) { pp_prefix = name; }

    // Set 'this' from multiblock container
    void SetMultiBlockPointer (MultiBlockContainer *mbc) { m_mbc = mbc; }

    // Public data copy for MB
    std::vector<amrex::Box> domain_p;
    MultiBlockContainer *m_mbc = nullptr;
    amrex::Vector<amrex::Vector<amrex::MultiFab> > vars_new;
    amrex::Vector<amrex::Vector<amrex::MultiFab> > vars_old;

    // Velocity time averaged field
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vel_t_avg;
    amrex::Vector<int> t_avg_cnt;

#endif
    std::string pp_prefix {"erf"};

    void fill_from_bndryregs (const amrex::Vector<amrex::MultiFab*>& mfs,
                              amrex::Real time);

#ifdef ERF_USE_NETCDF
    void fill_from_realbdy (const amrex::Vector<amrex::MultiFab*>& mfs,
                            amrex::Real time,
                            bool cons_only,
                            int icomp_cons,
                            int ncomp_cons);
#endif

#ifdef ERF_USE_NETCDF
    void init_from_wrfinput (int lev);
    void init_from_metgrid (int lev);
#endif // ERF_USE_NETCDF

#ifdef ERF_USE_WINDFARM
    void init_windfarm(int lev);
#endif

#ifdef ERF_USE_EB
    void MakeEBGeometry ();
    void make_eb_box ();
    void make_eb_cylinder ();
    void make_eb_regular ();
#endif

    // more flexible version of AverageDown() that lets you average down across multiple levels
    void AverageDownTo (int crse_lev, int scomp, int ncomp); // NOLINT

private:

    ///////////////////////////
    // private member functions
    ///////////////////////////

    // read in some parameters from inputs file
    void ReadParameters ();

    // set covered coarse cells/faces to be the average of overlying fine cells/faces
    void AverageDown ();

    void update_diffusive_arrays (int lev, const amrex::BoxArray& ba, const amrex::DistributionMapping& dm);

    void update_terrain_arrays (int lev, amrex::Real time);

    void Construct_ERFFillPatchers (int lev);

    void Define_ERFFillPatchers (int lev);

    void init1DArrays ();

    void init_bcs ();

    void init_custom (int lev);

    void init_uniform (int lev);

    void init_stuff (int lev, const amrex::BoxArray& ba, const amrex::DistributionMapping& dm,
                     amrex::Vector<amrex::MultiFab>& lev_new,
                     amrex::Vector<amrex::MultiFab>& lev_old);

    // Initialize the integrator object
    void initialize_integrator (int lev, amrex::MultiFab& cons_mf, amrex::MultiFab& vel_mf);

    // Initialize the bc routines
    void initialize_bcs        (int lev);

    // Initialize the microphysics object
    void initializeMicrophysics (const int&);

    // Compute a vector of new MultiFabs by copying from valid region and filling ghost cells
    //
    // NOTE: FillPatch takes in an empty MF, and returns cell-centered + velocities (not momenta)
    //
    void FillPatch (int lev, amrex::Real time,
                    const amrex::Vector<amrex::MultiFab*>& mfs_vel,
                    const amrex::Vector<amrex::MultiFab*>& mfs_mom,
                    bool fillset=true, bool cons_only=false);

    // Compute a new MultiFab by copying from valid region and filling ghost cells -
    void FillPatchMoistVars (int lev, amrex::MultiFab& mf);

    // Compute new multifabs by copying data from valid region and filling ghost cells.
    // Unlike FillPatch, FillIntermediatePatch will use the supplied multifabs instead of fine level data.
    // This is to support filling boundary cells at an intermediate time between old/new times
    // on the fine level when valid data at a specific time is already available (such as
    // at each RK stage when integrating between initial and final times at a given level).
    //
    // NOTE: FillIntermediatePatch takes in updated momenta, and returns both updated velocity and momenta
    //
    void FillIntermediatePatch (int lev, amrex::Real time,
                                const amrex::Vector<amrex::MultiFab*>& mfs_vel,
                                const amrex::Vector<amrex::MultiFab*>& mfs_mom,
                                int ng_cons, int ng_vel, bool cons_only, int icomp_cons, int ncomp_cons,
                                bool allow_most_bcs = true);

    // Fill all multifabs (and all components) in a vector of multifabs corresponding to the
    // grid variables defined in vars_old and vars_new just as FillCoarsePatch.
    void FillCoarsePatch (int lev, amrex::Real time);

    // advance a level by dt
    // includes a recursive call for finer levels
    void timeStep (int lev, amrex::Real time, int iteration);

    // advance a single level for a single time step
    void Advance (int lev, amrex::Real time, amrex::Real dt_lev, int iteration, int ncycle);

    //! Initialize HSE
    void initHSE ();
    void initHSE (int lev);

    //! Initialize Rayleigh damping profiles
    void initRayleigh ();

    //! Initialize sponge profiles
    void initSponge ();

    //! Set Rayleigh mean profiles from input sounding
    void setRayleighRefFromSounding (bool restarting);

    //! Set sponge mean profiles from input sounding
    void setSpongeRefFromSounding (bool restarting);

    // a wrapper for estTimeStep()
    void ComputeDt ();

    // get plotfile name
    [[nodiscard]] std::string PlotFileName (int lev) const;

    // set plotfile variables names
    static amrex::Vector<std::string> PlotFileVarNames (amrex::Vector<std::string> plot_var_names) ;

    // set which variables and derived quantities go into plotfiles
    void setPlotVariables (const std::string& pp_plot_var_names, amrex::Vector<std::string>& plot_var_names);
    // append variables to plot
    void appendPlotVariables (const std::string& pp_plot_var_names, amrex::Vector<std::string>& plot_var_names);

#ifdef ERF_USE_NETCDF
    //! Write plotfile using NETCDF
    void writeNCPlotFile (int lev, int which, const std::string& dir,
                          const amrex::Vector<const amrex::MultiFab*> &mf,
                          const amrex::Vector<std::string> &plot_var_names,
                          const amrex::Vector<int>& level_steps, amrex::Real time) const;

    //! Write checkpointFile using NetCdf
    void WriteNCCheckpointFile () const;

    //! Read checkpointFile for restart
    void ReadNCCheckpointFile ();

    //! Write MultiFab in NetCDF format
    static void WriteNCMultiFab (const amrex::FabArray<amrex::FArrayBox> &fab,
                                 const std::string& name,
                                 bool set_ghost = false) ;

    //! Read MultiFab in NetCDF format
    void ReadNCMultiFab (amrex::FabArray<amrex::FArrayBox> &mf,
                         const std::string &name,
                         int coordinatorProc = amrex::ParallelDescriptor::IOProcessorNumber(),
                         int allow_empty_mf = 0);

    //! Create 1D vertical column output for coupling
    void createNCColumnFile (int lev,
                             const std::string& colfile_name, amrex::Real xloc, amrex::Real yloc);

    // Copy from the NC*fabs into the MultiFabs holding the boundary data
    void init_from_wrfbdy (amrex::Vector<amrex::FArrayBox*> x_vel_lateral,
                           amrex::Vector<amrex::FArrayBox*> y_vel_lateral,
                           amrex::Vector<amrex::FArrayBox*> z_vel_lateral,
                           amrex::Vector<amrex::FArrayBox*> T_lateral);

    amrex::Real start_bdy_time;

    // *** *** FArrayBox's for holding the SURFACE data
    // amrex::IArrayBox NC_IVGTYP_fab; // Vegetation type (IVGTYP); Discrete numbers;
    // amrex::FArrayBox NC_z0_fab;     // Surface Roughness, z0 = z0 (IVGTYP)
    // amrex::FArrayBox NC_PSFC_fab;   // Surface pressure

    // TODO: Clarify the relation between SST and TSK
    // amrex::FArrayBox NC_SST_fab;    // Sea Surface Temperature; Defined even for land area
    // amrex::FArrayBox NC_TSK_fab;    // Surface Skin Temperature; Appears to be same as SST...

    // Vectors (over time) of Vector (over variables) of FArrayBoxs for holding the data read from the wrfbdy NetCDF file
    amrex::Vector<amrex::Vector<amrex::FArrayBox>> bdy_data_xlo;
    amrex::Vector<amrex::Vector<amrex::FArrayBox>> bdy_data_xhi;
    amrex::Vector<amrex::Vector<amrex::FArrayBox>> bdy_data_ylo;
    amrex::Vector<amrex::Vector<amrex::FArrayBox>> bdy_data_yhi;

    amrex::Real bdy_time_interval;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> lat_m, lon_m;
    amrex::Real Latitude;
    amrex::Real Longitude;
#endif // ERF_USE_NETCDF

    // Struct for working with the sounding data we take as an input
    InputSoundingData input_sounding_data;

    // Struct for working with the sponge data we take as an input
    InputSpongeData input_sponge_data;

    // write checkpoint file to disk
    void WriteCheckpointFile () const;

    // read checkpoint file from disk
    void ReadCheckpointFile ();

    // Read the file passed to amr.restart and use it as an initial condition for
    // the current simulation. Supports a different number of components and
    // ghost cells.
    void InitializeFromFile ();

    // Initialize the new-time data at a level from the initial_data MultiFab
    void InitializeLevelFromData (int lev, const amrex::MultiFab& initial_data);

    // utility to skip to next line in Header
    static void GotoNextLine (std::istream& is);

    // Single level functions called by advance()
    void post_update (amrex::MultiFab& state_mf, amrex::Real time, const amrex::Geometry& geom);
    void fill_rhs (amrex::MultiFab& rhs_mf, const amrex::MultiFab& state_mf, amrex::Real time, const amrex::Geometry& geom);

    ////////////////
    // private data members

    // These are the minimum and maximum values of k in valid boxes *at this level*
    amrex::Vector<int> min_k_at_level;
    amrex::Vector<int> max_k_at_level;

    std::unique_ptr<ProblemBase> prob = nullptr;  // problem-specific functions

    amrex::Vector<int> num_boxes_at_level;                   // how many boxes specified at each level by tagging criteria
    amrex::Vector<int> num_files_at_level;                   // how many wrfinput files specified at each level
    amrex::Vector<amrex::Vector<amrex::Box>> boxes_at_level; //      the boxes specified at each level by tagging criteria

    amrex::Vector<int> istep;      // which step?
    amrex::Vector<int> nsubsteps;  // how many substeps on each level?

    // keep track of old time, new time, and time step at each level
    amrex::Vector<amrex::Real> t_new;
    amrex::Vector<amrex::Real> t_old;
    amrex::Vector<amrex::Real> dt;
    amrex::Vector<long> dt_mri_ratio;

    // array of multifabs to store the solution at each level of refinement
    // after advancing a level we use "swap".
#ifndef ERF_USE_MULTIBLOCK
    amrex::Vector<amrex::Vector<amrex::MultiFab> > vars_new;
    amrex::Vector<amrex::Vector<amrex::MultiFab> > vars_old;

    // Velocity time averaged field
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vel_t_avg;
    amrex::Vector<amrex::Real> t_avg_cnt;
#endif
    amrex::Vector<std::unique_ptr<MRISplitIntegrator<amrex::Vector<amrex::MultiFab> > > > mri_integrator_mem;

    // Vector over levels of routines to impose physical boundary conditions
    amrex::Vector<std::unique_ptr<ERFPhysBCFunct_cons>> physbcs_cons;
    amrex::Vector<std::unique_ptr<ERFPhysBCFunct_u>> physbcs_u;
    amrex::Vector<std::unique_ptr<ERFPhysBCFunct_v>> physbcs_v;
    amrex::Vector<std::unique_ptr<ERFPhysBCFunct_w>> physbcs_w;
    amrex::Vector<std::unique_ptr<ERFPhysBCFunct_w_no_terrain>> physbcs_w_no_terrain;

    // Store primitive variable for MOST BC
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> Theta_prim;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> Qv_prim;

    // Scratch space for time integrator
    amrex::Vector<amrex::MultiFab> rU_old;
    amrex::Vector<amrex::MultiFab> rU_new;
    amrex::Vector<amrex::MultiFab> rV_old;
    amrex::Vector<amrex::MultiFab> rV_new;
    amrex::Vector<amrex::MultiFab> rW_old;
    amrex::Vector<amrex::MultiFab> rW_new;

    std::unique_ptr<Microphysics> micro;
    amrex::Vector<amrex::Vector<amrex::MultiFab*>> qmoist; // (lev,ncomp) This has up to 8 components: qt, qv, qc, qi, qp, qr, qs, qg

    amrex::Vector<amrex::MultiFab> Nturb;
    amrex::Vector<amrex::MultiFab> vars_fitch; // Vabs, Vabsdt, dudt, dvdt, dTKEdt
    amrex::Vector<amrex::MultiFab> vars_ewp; // dudt, dvdt, dTKEdt
    amrex::Real hub_height, rotor_dia, thrust_coeff_standing, nominal_power;
    amrex::Vector<amrex::Real> wind_speed, thrust_coeff, power;


    LandSurface lsm;
    amrex::Vector<amrex::Vector<amrex::MultiFab*>> lsm_data; // (lev,ncomp) Components: theta, q1, q2
    amrex::Vector<amrex::Vector<amrex::MultiFab*>> lsm_flux; // (lev,ncomp) Components: theta, q1, q2

#if defined(ERF_USE_RRTMGP)
    Radiation rad;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> qheating_rates;  // radiation heating rate source terms

    // Containers for additional SLM inputs
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> sw_lw_fluxes; // Direct SW (visible, NIR), Diffuse SW (visible, NIR), LW flux
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> solar_zenith; // Solar zenith angle
#endif

    // Fillpatcher classes for coarse-fine boundaries
    int cf_width{0};
    int cf_set_width{0};
    amrex::Vector<ERFFillPatcher> FPr_c;
    amrex::Vector<ERFFillPatcher> FPr_u;
    amrex::Vector<ERFFillPatcher> FPr_v;
    amrex::Vector<ERFFillPatcher> FPr_w;

    // Diffusive stresses and Smag
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> Tau11_lev, Tau22_lev, Tau33_lev;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> Tau12_lev, Tau21_lev;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> Tau13_lev, Tau31_lev;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> Tau23_lev, Tau32_lev;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> eddyDiffs_lev;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> SmnSmn_lev;

    // Sea Surface Temps and Land Masks (lev, ntimes)
    amrex::Vector<amrex::Vector<std::unique_ptr<amrex::MultiFab>>>  sst_lev;
    amrex::Vector<amrex::Vector<std::unique_ptr<amrex::iMultiFab>>> lmask_lev;

    // Other SFS terms
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> SFS_hfx1_lev, SFS_hfx2_lev, SFS_hfx3_lev;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> SFS_diss_lev;

    // Terrain / grid stretching
    amrex::Vector<amrex::Real> zlevels_stag; // nominal height levels
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> z_phys_nd;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> z_phys_cc;

    amrex::Vector<std::unique_ptr<amrex::MultiFab>>   detJ_cc;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>>   ax;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>>   ay;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>>   az;

    amrex::Vector<std::unique_ptr<amrex::MultiFab>> z_phys_nd_src;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>>   detJ_cc_src;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>>   ax_src;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>>   ay_src;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>>   az_src;

    amrex::Vector<std::unique_ptr<amrex::MultiFab>> z_phys_nd_new;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>>   detJ_cc_new;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>>   ax_new;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>>   ay_new;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>>   az_new;

    amrex::Vector<std::unique_ptr<amrex::MultiFab>> z_t_rk;

    amrex::Vector<std::unique_ptr<amrex::MultiFab>> mapfac_m;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> mapfac_u;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> mapfac_v;

    amrex::Vector<amrex::MultiFab> base_state;
    amrex::Vector<amrex::MultiFab> base_state_new;

    // array of flux registers
    amrex::Vector<amrex::YAFluxRegister*> advflux_reg;

    // A BCRec is essentially a 2*DIM integer array storing the boundary
    // condition type at each lo/hi walls in each direction. We have one BCRec
    // for each component of the cell-centered variables and each velocity component.
    amrex::Vector           <amrex::BCRec> domain_bcs_type;
    amrex::Gpu::DeviceVector<amrex::BCRec> domain_bcs_type_d;

    // We store these so that we can print them out in the job_info file
    amrex::Array<std::string,2*AMREX_SPACEDIM> domain_bc_type;

    // These hold the Dirichlet values at walls which need them ...
    amrex::Array<amrex::Array<amrex::Real, AMREX_SPACEDIM*2>, AMREX_SPACEDIM+NVAR_max> m_bc_extdir_vals;

    // These hold the Neumann values at walls which need them ...
    amrex::Array<amrex::Array<amrex::Real, AMREX_SPACEDIM*2>, AMREX_SPACEDIM+NVAR_max> m_bc_neumann_vals;

    // These are the "physical" boundary condition types (e.g. "inflow")
    amrex::GpuArray<ERF_BC, AMREX_SPACEDIM*2> phys_bc_type;

    // These are the masks for the thin immersed body
    amrex::Vector<std::unique_ptr<amrex::iMultiFab>> xflux_imask;
    amrex::Vector<std::unique_ptr<amrex::iMultiFab>> yflux_imask;
    amrex::Vector<std::unique_ptr<amrex::iMultiFab>> zflux_imask;
    //amrex::Vector<std::unique_ptr<amrex::iMultiFab>> overset_imask;

    // These are the body forces that result from the thin immersed body
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> thin_xforce;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> thin_yforce;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> thin_zforce;

    int last_plot_file_step_1;
    int last_plot_file_step_2;

    int last_check_file_step;
    int plot_file_on_restart = 1;

    ////////////////
    // runtime parameters

    // maximum number of steps and stop time
    int max_step = std::numeric_limits<int>::max();
    amrex::Real start_time = 0.0;
    amrex::Real  stop_time = std::numeric_limits<amrex::Real>::max();

    // if >= 0 we restart from a checkpoint
    std::string restart_chkfile = "";

    // Time step controls
    static amrex::Real cfl;
    static amrex::Real init_shrink;
    static amrex::Real change_max;

    // Fixed dt for level 0 timesteps (only used if positive)
    static amrex::Real fixed_dt;
    static amrex::Real fixed_fast_dt;
    static int fixed_mri_dt_ratio;

    // how often each level regrids the higher levels of refinement
    // (after a level advances that many time steps)
    int regrid_int = -1;

    // plotfile prefix and frequency
    std::string plot_file_1 {"plt_1_"};
    std::string plot_file_2 {"plt_2_"};
    int m_plot_int_1 = -1;
    int m_plot_int_2 = -1;
    amrex::Real m_plot_per_1 = -1.0;
    amrex::Real m_plot_per_2 = -1.0;

    bool plot_lsm = false;

    // other sampling output control
    int profile_int = -1;
    bool cc_profiles = true;

    // Checkpoint type, prefix and frequency
    std::string check_file {"chk"};
    std::string check_type {"native"};
    std::string restart_type {"native"};
    int m_check_int = -1;
    amrex::Real m_check_per = -1.0;

    amrex::Vector<std::string> plot_var_names_1;
    amrex::Vector<std::string> plot_var_names_2;
    const amrex::Vector<std::string> cons_names     {"density", "rhotheta", "rhoKE", "rhoQKE", "rhoadv_0",
                                                     "rhoQ1", "rhoQ2", "rhoQ3",
                                                     "rhoQ4", "rhoQ5", "rhoQ6"};

    // Note that the order of variable names here must match the order in Derive.cpp
    const amrex::Vector<std::string> derived_names {"soundspeed", "temp", "theta", "KE", "QKE", "scalar",
                                                    "vorticity_x","vorticity_y","vorticity_z",
                                                    "magvel", "divU",
                                                    "pres_hse", "dens_hse", "pressure", "pert_pres", "pert_dens", "eq_pot_temp", "num_turb",
                                                    "dpdx", "dpdy", "pres_hse_x", "pres_hse_y",
                                                    "z_phys", "detJ" , "mapfac", "lat_m", "lon_m",
                                                    // Time averaged velocity
                                                    "u_t_avg", "v_t_avg", "w_t_avg", "umag_t_avg",
                                                    // eddy viscosity
                                                    "Kmv","Kmh",
                                                    // eddy diffusivity of heat
                                                    "Khv","Khh",
                                                    // mynn pbl lengthscale
                                                    "Lpbl",
                                                    // moisture vars
                                                    "qt", "qv", "qc", "qi", "qp", "qrain", "qsnow", "qgraup", "rain_accum", "snow_accum", "graup_accum"
#ifdef ERF_COMPUTE_ERROR
                                                    ,"xvel_err", "yvel_err", "zvel_err", "pp_err"
#endif
                                                   };

    // algorithm choices
    static SolverChoice solverChoice;

#ifdef ERF_USE_PARTICLES
    // Particle container with all particle species
    ParticleData particleData;

    // variables and functions for tracers particles
    bool  m_use_tracer_particles; /*!< tracer particles that advect with flow */
    bool  m_use_hydro_particles;  /*!< tracer particles that fall with gravity */

    /*! Read tracer and hydro particles parameters */
    void readTracersParams();

    /*! Initialize tracer and hydro particles */
    void initializeTracers ( amrex::ParGDBBase*,
                             const amrex::Vector<std::unique_ptr<amrex::MultiFab>>& );
    /*! Evolve tracers and hydro particles */
    void evolveTracers( int,
                        amrex::Real,
                        amrex::Vector<amrex::Vector<amrex::MultiFab>>&,
                        const amrex::Vector<std::unique_ptr<amrex::MultiFab>>& );

#endif

    static int verbose;

    // Diagnostic output interval
    static int sum_interval;
    static amrex::Real sum_per;

    // Native or NetCDF
    static std::string plotfile_type;

    // init_type:  "ideal", "real", "input_sounding", "metgrid" or ""
    static std::string init_type;

    // sponge_type:  "input_sponge"
    static std::string sponge_type;

    // use_real_bcs: only true if 1) ( (init_type == real) or (init_type == metgrid) )
    //                        AND 2) we want to use the bc's from the WRF bdy file
    static bool use_real_bcs;

    // NetCDF initialization (wrfinput) file
    static amrex::Vector<amrex::Vector<std::string>> nc_init_file;

    // NetCDF initialization (wrfbdy/met_em) file
    static std::string nc_bdy_file;
    int real_width{0};
    int real_set_width{0};

    // Text input_sounding file
    static std::string input_sounding_file;

    // Text input_sponge file
    static std::string input_sponge_file;

    // Flag to trigger initialization from input_sounding like WRF's ideal.exe
    // used with init_type == "input_sounding"
    static bool init_sounding_ideal;

    // 1D CDF output (for ingestion in AMR-Wind)
    static int         output_1d_column;
    static int         column_interval;
    static amrex::Real column_per;
    static amrex::Real column_loc_x;
    static amrex::Real column_loc_y;
    static std::string column_file_name;

    // 2D BndryRegister output (for ingestion in AMR-Wind)
    static int         output_bndry_planes;
    static int         bndry_output_planes_interval;
    static amrex::Real bndry_output_planes_per;
    static amrex::Real bndry_output_planes_start_time;

    // 2D BndryRegister input
    static int          input_bndry_planes;

    static int ng_dens_hse;
    static int ng_pres_hse;

    // Custom source terms
    amrex::Vector<           amrex::Vector<amrex::Real> > h_rhotheta_src;
    amrex::Vector<amrex::Gpu::DeviceVector<amrex::Real> > d_rhotheta_src;

    amrex::Vector<           amrex::Vector<amrex::Real> > h_rhoqt_src;
    amrex::Vector<amrex::Gpu::DeviceVector<amrex::Real> > d_rhoqt_src;

    amrex::Vector<           amrex::Vector<amrex::Real> > h_w_subsid;
    amrex::Vector<amrex::Gpu::DeviceVector<amrex::Real> > d_w_subsid;

    amrex::Vector<           amrex::Vector<amrex::Real> > h_u_geos;
    amrex::Vector<amrex::Gpu::DeviceVector<amrex::Real> > d_u_geos;

    amrex::Vector<           amrex::Vector<amrex::Real> > h_v_geos;
    amrex::Vector<amrex::Gpu::DeviceVector<amrex::Real> > d_v_geos;

    // This is a vector over levels of vectors across quantities of Vectors
    amrex::Vector<amrex::Vector<           amrex::Vector<amrex::Real> > > h_rayleigh_ptrs;

    // This is a vector over levels of vectors across quantities of Vectors
    amrex::Vector<amrex::Vector<           amrex::Vector<amrex::Real> > > h_sponge_ptrs;

    // This is a vector over levels of vectors across quantities of DeviceVectors
    amrex::Vector<amrex::Vector<amrex::Gpu::DeviceVector<amrex::Real> > > d_rayleigh_ptrs;

    // This is a vector over levels of vectors across quantities of DeviceVectors
    amrex::Vector<amrex::Vector<amrex::Gpu::DeviceVector<amrex::Real> > > d_sponge_ptrs;

    amrex::Vector<amrex::Real> h_havg_density;
    amrex::Vector<amrex::Real> h_havg_temperature;
    amrex::Vector<amrex::Real> h_havg_pressure;
    amrex::Vector<amrex::Real> h_havg_qv;
    amrex::Vector<amrex::Real> h_havg_qc;

    amrex::Gpu::DeviceVector<amrex::Real> d_havg_density;
    amrex::Gpu::DeviceVector<amrex::Real> d_havg_temperature;
    amrex::Gpu::DeviceVector<amrex::Real> d_havg_pressure;
    amrex::Gpu::DeviceVector<amrex::Real> d_havg_qv;
    amrex::Gpu::DeviceVector<amrex::Real> d_havg_qc;

    void refinement_criteria_setup ();

    std::unique_ptr<WriteBndryPlanes> m_w2d  = nullptr;
    std::unique_ptr<ReadBndryPlanes>  m_r2d  = nullptr;
    std::unique_ptr<ABLMost>          m_most = nullptr;

    //
    // Holds info for dynamically generated tagging criteria
    //
    static amrex::Vector<amrex::AMRErrorTag> ref_tags;

    //
    // Build a mask that zeroes out values on a coarse level underlying
    //     grids on the next finest level
    //
    amrex::MultiFab fine_mask;

    amrex::Real dz_min;

    static AMREX_FORCE_INLINE
    int
    ComputeGhostCells (const AdvChoice& advChoice, bool use_num_diff)
    {
        if (use_num_diff)
        {
            return 3;
        } else {
            if (
                 (advChoice.dycore_horiz_adv_type    == AdvType::Centered_6th)
              || (advChoice.dycore_vert_adv_type     == AdvType::Centered_6th)
              || (advChoice.moistscal_horiz_adv_type == AdvType::Centered_6th)
              || (advChoice.moistscal_horiz_adv_type == AdvType::Centered_6th)
              || (advChoice.dryscal_horiz_adv_type   == AdvType::Centered_6th)
              || (advChoice.dryscal_vert_adv_type    == AdvType::Centered_6th) )
            { return 3; }
            else if (
                 (advChoice.dycore_horiz_adv_type    == AdvType::Upwind_5th)
              || (advChoice.dycore_vert_adv_type     == AdvType::Upwind_5th)
              || (advChoice.moistscal_horiz_adv_type == AdvType::Upwind_5th)
              || (advChoice.moistscal_horiz_adv_type == AdvType::Upwind_5th)
              || (advChoice.dryscal_horiz_adv_type   == AdvType::Upwind_5th)
              || (advChoice.dryscal_vert_adv_type    == AdvType::Upwind_5th) )
            { return 3; }
            else if (
                  (advChoice.dryscal_horiz_adv_type   == AdvType::Weno_5)
               || (advChoice.dryscal_vert_adv_type    == AdvType::Weno_5)
               || (advChoice.moistscal_horiz_adv_type == AdvType::Weno_5)
               || (advChoice.moistscal_vert_adv_type  == AdvType::Weno_5)
               || (advChoice.moistscal_horiz_adv_type == AdvType::Weno_5Z)
               || (advChoice.moistscal_vert_adv_type  == AdvType::Weno_5Z)
               || (advChoice.dryscal_horiz_adv_type   == AdvType::Weno_5Z)
               || (advChoice.dryscal_vert_adv_type    == AdvType::Weno_5Z) )
            { return 3; }
            else
            { return 2; }
        }
    }

    AMREX_FORCE_INLINE
    amrex::YAFluxRegister* getAdvFluxReg (int lev)
    {
        return advflux_reg[lev];
    }

    AMREX_FORCE_INLINE
    std::ostream&
    DataLog (int i)
    {
        return *datalog[i];
    }

    AMREX_FORCE_INLINE
    int
    NumDataLogs () noexcept
    {
        return datalog.size();
    }

    AMREX_FORCE_INLINE
    std::ostream&
    SamplePointLog (int i)
    {
        return *sampleptlog[i];
    }

    AMREX_FORCE_INLINE
    int
    NumSamplePointLogs () noexcept
    {
        return sampleptlog.size();
    }

    AMREX_FORCE_INLINE
    std::ostream&
    SampleLineLog (int i)
    {
        return *samplelinelog[i];
    }

    AMREX_FORCE_INLINE
    int
    NumSampleLineLogs () noexcept
    {
        return samplelinelog.size();
    }

    amrex::IntVect&
    SamplePoint (int i)
    {
        return samplepoint[i];
    }

    AMREX_FORCE_INLINE
    int
    NumSamplePoints () noexcept
    {
        return samplepoint.size();
    }

    amrex::IntVect&
    SampleLine (int i)
    {
        return sampleline[i];
    }

    AMREX_FORCE_INLINE
    int
    NumSampleLines () noexcept
    {
        return sampleline.size();
    }

    static amrex::Real startCPUTime;
    static amrex::Real previousCPUTimeUsed;

    static amrex::Real
    getCPUTime ()
    {
      int numCores = amrex::ParallelDescriptor::NProcs();
#ifdef _OPENMP
      numCores = numCores * omp_get_max_threads();
#endif

      amrex::Real T =
        numCores * (amrex::ParallelDescriptor::second() - startCPUTime) +
        previousCPUTimeUsed;

      return T;
    }

    void setRecordDataInfo (int i, const std::string& filename) // NOLINT
    {
        if (amrex::ParallelDescriptor::IOProcessor())
        {
            datalog[i] = std::make_unique<std::fstream>();
            datalog[i]->open(filename.c_str(),std::ios::out|std::ios::app);
            if (!datalog[i]->good()) {
                amrex::FileOpenFailed(filename);
            }
        }
        amrex::ParallelDescriptor::Barrier("ERF::setRecordDataInfo");
    }

    void setRecordSamplePointInfo (int i, int lev, amrex::IntVect& cell, const std::string& filename) // NOLINT
    {
        amrex::MultiFab dummy(grids[lev],dmap[lev],1,0);
        for (amrex::MFIter mfi(dummy); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            if (bx.contains(cell)) {
                sampleptlog[i] = std::make_unique<std::fstream>();
                sampleptlog[i]->open(filename.c_str(),std::ios::out|std::ios::app);
                if (!sampleptlog[i]->good()) {
                    amrex::FileOpenFailed(filename);
                }
            }
        }
        amrex::ParallelDescriptor::Barrier("ERF::setRecordSamplePointInfo");
    }

    void setRecordSampleLineInfo (int i, int lev, amrex::IntVect& cell, const std::string& filename) // NOLINT
    {
        amrex::MultiFab dummy(grids[lev],dmap[lev],1,0);
        for (amrex::MFIter mfi(dummy); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            if (bx.contains(cell)) {
                samplelinelog[i] = std::make_unique<std::fstream>();
                samplelinelog[i]->open(filename.c_str(),std::ios::out|std::ios::app);
                if (!samplelinelog[i]->good()) {
                    amrex::FileOpenFailed(filename);
                }
            }
        }
        amrex::ParallelDescriptor::Barrier("ERF::setRecordSampleLineInfo");
    }

    amrex::Vector<std::unique_ptr<std::fstream> > datalog;
    amrex::Vector<std::string> datalogname;

    amrex::Vector<std::unique_ptr<std::fstream> > sampleptlog;
    amrex::Vector<std::string> sampleptlogname;
    amrex::Vector<amrex::IntVect> samplepoint;

    amrex::Vector<std::unique_ptr<std::fstream> > samplelinelog;
    amrex::Vector<std::string> samplelinelogname;
    amrex::Vector<amrex::IntVect> sampleline;

    //! The filename of the ith datalog file.
    [[nodiscard]] std::string DataLogName (int i) const noexcept { return datalogname[i]; }

    //! The filename of the ith sampleptlog file.
    [[nodiscard]] std::string SamplePointLogName (int i) const noexcept { return sampleptlogname[i]; }

    //! The filename of the ith samplelinelog file.
    [[nodiscard]] std::string SampleLineLogName (int i) const noexcept { return samplelinelogname[i]; }

    amrex::Vector<std::unique_ptr<amrex::FabFactory<amrex::FArrayBox> > > m_factory;

    [[nodiscard]] amrex::FabFactory<amrex::FArrayBox> const&
    Factory (int lev) const noexcept { return *m_factory[lev]; }
#ifdef ERF_USE_EB
    [[nodiscard]] amrex::EBFArrayBoxFactory const&
    EBFactory (int lev) const noexcept {
        return static_cast<amrex::EBFArrayBoxFactory const&>(*m_factory[lev]);
    }

    [[nodiscard]] static int nghost_eb_basic ()
    { return 5; }

    // We need 5 for doing StateRedistribution; otherwise 4 would be enough
    [[nodiscard]] static int nghost_eb_volume ()
    { return 5; }

    [[nodiscard]] static int nghost_eb_full ()
    { return 4; }
#endif

public:
    void writeJobInfo (const std::string& dir) const;
    static void writeBuildInfo (std::ostream& os);

    static void print_banner(MPI_Comm /*comm*/, std::ostream& /*out*/);
    static void print_usage(MPI_Comm /*comm*/, std::ostream& /*out*/);
    static void print_error(MPI_Comm /*comm*/, const std::string& msg);
    static void print_summary(std::ostream&);
    static void print_tpls(std::ostream& /*out*/);
};

#endif
