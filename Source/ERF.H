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
#include <AMReX_Print.H>

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_YAFluxRegister.H>
#include <AMReX_ErrorList.H>
#include <AMReX_EBFabFactory.H>
#include <AMReX_EBMultiFabUtil.H>

#include <ERF_EB.H>

#ifdef ERF_USE_FFT
#include <AMReX_FFT_Poisson.H>
#endif

#ifdef AMREX_MEM_PROFILING
#include <AMReX_MemProfiler.H>
#endif

#include "ERF_ProbCommon.H"

#include <ERF_IndexDefines.H>
#include <ERF_DataStruct.H>
#include <ERF_TurbPertStruct.H>
#include <ERF_InputSoundingData.H>
#include <ERF_InputSpongeData.H>
#include <ERF_SurfaceLayer.H>
#include <ERF_Derive.H>
#include <ERF_ReadBndryPlanes.H>
#include <ERF_WriteBndryPlanes.H>
#include <ERF_MRI.H>
#include <ERF_PhysBCFunct.H>
#include <ERF_FillPatcher.H>
#include <ERF_SampleData.H>
#include <ERF_ForestDrag.H>

#ifdef ERF_USE_PARTICLES
#include "ERF_ParticleData.H"
#endif

#include "ERF_EulerianMicrophysics.H"
#include "ERF_LagrangianMicrophysics.H"
#include "ERF_LandSurface.H"

#ifdef ERF_USE_WINDFARM
#include "ERF_WindFarm.H"
#endif

#include <ERF_RadiationModels.H>

#ifdef ERF_USE_SHOC
#include "ERF_ShocInterface.H"
#endif

#ifdef ERF_USE_P3
#include "ERF_P3Interface.H"
#endif

#include <iostream>

#ifdef AMREX_LAZY
#include <AMReX_Lazy.H>
#endif

#ifndef AMREX_USE_MPI
using amrex::MPI_COMM_WORLD;
using amrex::MPI_Comm;
#endif

#ifdef ERF_USE_MULTIBLOCK
class MultiBlockContainer;
#endif

/**
 * Enum of possible interpolation types between coarse/fine
*/
AMREX_ENUM(StateInterpType,
    FullState, Perturbational
);

/**
 * Enum of possible plotfile types
*/
AMREX_ENUM(PlotFileType,
    None, Amrex, Netcdf
);

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

    void ERF_shared ();

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

    // Refinement criteria for hurricanes
    void HurricaneTracker (int lev,
                           amrex::Real time,
                           const amrex::MultiFab& cc_vel,
                           const amrex::Real velmag_threshold,
                           amrex::TagBoxArray* tags = nullptr);

    bool FindInitialEye (int lev,
                         const amrex::MultiFab& cc_vel,
                         const amrex::Real velmag_threshold,
                         amrex::Real& eye_x, amrex::Real& eye_y);

    // Hurricane track line output as VTKPolyline
    amrex::Vector<std::array<amrex::Real, 2>> hurricane_track_xy;
    amrex::Vector<std::array<amrex::Real, 2>> hurricane_eye_track_xy,
                                              hurricane_eye_track_latlon,
                                              hurricane_maxvel_vs_time,
                                              hurricane_tracker_circle;
    amrex::Vector<amrex::MultiFab> weather_forecast_data_1, weather_forecast_data_2;
    amrex::Vector<amrex::Vector<amrex::MultiFab>> forecast_state_1,
                                                  forecast_state_2,
                                                  forecast_state_interp;
    std::string MakeVTKFilename(int nstep);
    std::string MakeVTKFilename_TrackerCircle(int nstep);
    std::string MakeVTKFilename_EyeTracker_xy(int nstep);
    std::string MakeFilename_EyeTracker_latlon(int nstep);
    std::string MakeFilename_EyeTracker_maxvel(int nstep);
    void WriteVTKPolyline(const std::string& filename,
                          amrex::Vector<std::array<amrex::Real, 2>>& points_xy);

    void WriteLinePlot(const std::string& filename,
                       amrex::Vector<std::array<amrex::Real, 2>>& points_xy);

    // Initialize multilevel data
    void InitData ();

    // Initialize multilevel data before MultiBlock
    void InitData_pre ();

    // Initialize multilevel data after MultiBlock
    void InitData_post ();

    // Interpolate 2D arrays lat_m, lon_m, sinPhi, cosPhi
    void Interp2DArrays (int lev, const amrex::BoxArray& my_ba2d, const amrex::DistributionMapping& my_dm);

    void WriteMyEBSurface ();

    // Compute the divergence -- whether EB, no-terrain, flat terrain or general terrain
    void compute_divergence (int lev, amrex::MultiFab& rhs, amrex::Array<amrex::MultiFab const*,AMREX_SPACEDIM> rho0_u_const,
                             amrex::Geometry const& geom_at_lev);

    // Project the velocities to be divergence-free -- this is only relevant if anelastic == 1
    void project_initial_velocity (int lev, amrex::Real time, amrex::Real dt);
    void project_momenta (int lev, amrex::Real l_time, amrex::Real l_dt, amrex::Vector<amrex::MultiFab >& vars);

    // Project the velocities to be divergence-free with a thin body
    void project_velocity_tb (int lev, amrex::Real dt, amrex::Vector<amrex::MultiFab >& vars);

    // Calculate wall distance by solving a Poisson equation
    void poisson_wall_dist (int lev);

    void make_subdomains(const amrex::BoxList& ba, amrex::Vector<amrex::BoxArray>& bins);

#ifdef ERF_USE_FFT
    void solve_with_fft (int lev, const amrex::Box& subdomain,
                         amrex::MultiFab& rhs, amrex::MultiFab& p,
                         amrex::Array<amrex::MultiFab,AMREX_SPACEDIM>& fluxes);
#endif
    void solve_with_gmres (int lev, const amrex::Box& subdomain,
                           amrex::MultiFab& rhs, amrex::MultiFab& p,
                           amrex::Array<amrex::MultiFab,AMREX_SPACEDIM>& fluxes,
                           amrex::MultiFab& ax_sub, amrex::MultiFab& ay_sub,
                           amrex::MultiFab& az_sub,
                           amrex::MultiFab& /*dJ_sub*/, amrex::MultiFab& znd_sub);

    void ImposeBCsOnPhi (int lev, amrex::MultiFab& phi, const amrex::Box& subdomain);

    // Init (NOT restart or regrid)
    void init_only (int lev, amrex::Real time);

    // Restart
    void restart ();

    // Check for NaNs and other bad stuff
    void check_state_for_nans(amrex::MultiFab const& S);
    void check_vels_for_nans(amrex::MultiFab const& xvel, amrex::MultiFab const& yvel, amrex::MultiFab const& zvel);
    void check_for_negative_theta (amrex::MultiFab& S);
    void check_for_low_temp       (amrex::MultiFab& S);

    // Is it time to write a plotfile or checkpoint?
    bool writeNow (const amrex::Real cur_time, const int nstep, const int plot_int,
                   const amrex::Real plot_per, const amrex::Real dt_0, amrex::Real& last_file_time);

    // Called after every level 0 timestep
    void post_timestep (int nstep, amrex::Real time, amrex::Real dt_lev);

    // Diagnostics
    void sum_integrated_quantities (amrex::Real time);
    void sum_derived_quantities (amrex::Real time);
    void sum_energy_quantities (amrex::Real time);

    void write_1D_profiles (amrex::Real time);
    void write_1D_profiles_stag (amrex::Real time);

    amrex::Real cloud_fraction (amrex::Real time);

    // Fill the physical boundary conditions for cell-centered velocity (diagnostic only)
    void FillBdyCCVels (amrex::Vector<amrex::MultiFab>& mf_cc_vel, int levc=0);

    void sample_points (int lev, amrex::Real time, amrex::IntVect cell, amrex::MultiFab& mf);
    void sample_lines (int lev, amrex::Real time, amrex::IntVect cell, amrex::MultiFab& mf);

    void derive_diag_profiles (amrex::Real time,
                               amrex::Gpu::HostVector<amrex::Real>& h_avg_u   , amrex::Gpu::HostVector<amrex::Real>& h_avg_v,
                               amrex::Gpu::HostVector<amrex::Real>& h_avg_w   , amrex::Gpu::HostVector<amrex::Real>& h_avg_rho,
                               amrex::Gpu::HostVector<amrex::Real>& h_avg_th  , amrex::Gpu::HostVector<amrex::Real>& h_avg_ksgs,
                               amrex::Gpu::HostVector<amrex::Real>& h_avg_Kmv , amrex::Gpu::HostVector<amrex::Real>& h_avg_Khv,
                               amrex::Gpu::HostVector<amrex::Real>& h_avg_qv  , amrex::Gpu::HostVector<amrex::Real>& h_avg_qc,
                               amrex::Gpu::HostVector<amrex::Real>& h_avg_qr  ,
                               amrex::Gpu::HostVector<amrex::Real>& h_avg_wqv , amrex::Gpu::HostVector<amrex::Real>& h_avg_wqc,
                               amrex::Gpu::HostVector<amrex::Real>& h_avg_wqr ,
                               amrex::Gpu::HostVector<amrex::Real>& h_avg_qi  , amrex::Gpu::HostVector<amrex::Real>& h_avg_qs,
                               amrex::Gpu::HostVector<amrex::Real>& h_avg_qg  ,
                               amrex::Gpu::HostVector<amrex::Real>& h_avg_uu  , amrex::Gpu::HostVector<amrex::Real>& h_avg_uv,
                               amrex::Gpu::HostVector<amrex::Real>& h_avg_uw,
                               amrex::Gpu::HostVector<amrex::Real>& h_avg_vv  , amrex::Gpu::HostVector<amrex::Real>& h_avg_vw,
                               amrex::Gpu::HostVector<amrex::Real>& h_avg_ww,
                               amrex::Gpu::HostVector<amrex::Real>& h_avg_uth , amrex::Gpu::HostVector<amrex::Real>& h_avg_vth,
                               amrex::Gpu::HostVector<amrex::Real>& h_avg_wth, amrex::Gpu::HostVector<amrex::Real>& h_avg_thth,
                               amrex::Gpu::HostVector<amrex::Real>& h_avg_ku, amrex::Gpu::HostVector<amrex::Real>& h_avg_kv,
                               amrex::Gpu::HostVector<amrex::Real>& h_avg_kw,
                               amrex::Gpu::HostVector<amrex::Real>& h_avg_p,
                               amrex::Gpu::HostVector<amrex::Real>& h_avg_pu, amrex::Gpu::HostVector<amrex::Real>& h_avg_pv,
                               amrex::Gpu::HostVector<amrex::Real>& h_avg_pw, amrex::Gpu::HostVector<amrex::Real>& h_avg_wthv);
    void derive_diag_profiles_stag (amrex::Real time,
                                    amrex::Gpu::HostVector<amrex::Real>& h_avg_u   , amrex::Gpu::HostVector<amrex::Real>& h_avg_v,
                                    amrex::Gpu::HostVector<amrex::Real>& h_avg_w   , amrex::Gpu::HostVector<amrex::Real>& h_avg_rho,
                                    amrex::Gpu::HostVector<amrex::Real>& h_avg_th  , amrex::Gpu::HostVector<amrex::Real>& h_avg_ksgs,
                                    amrex::Gpu::HostVector<amrex::Real>& h_avg_Kmv , amrex::Gpu::HostVector<amrex::Real>& h_avg_Khv,
                                    amrex::Gpu::HostVector<amrex::Real>& h_avg_qv  , amrex::Gpu::HostVector<amrex::Real>& h_avg_qc,
                                    amrex::Gpu::HostVector<amrex::Real>& h_avg_qr  ,
                                    amrex::Gpu::HostVector<amrex::Real>& h_avg_wqv , amrex::Gpu::HostVector<amrex::Real>& h_avg_wqc,
                                    amrex::Gpu::HostVector<amrex::Real>& h_avg_wqr ,
                                    amrex::Gpu::HostVector<amrex::Real>& h_avg_qi  , amrex::Gpu::HostVector<amrex::Real>& h_avg_qs,
                                    amrex::Gpu::HostVector<amrex::Real>& h_avg_qg  ,
                                    amrex::Gpu::HostVector<amrex::Real>& h_avg_uu  , amrex::Gpu::HostVector<amrex::Real>& h_avg_uv,
                                    amrex::Gpu::HostVector<amrex::Real>& h_avg_uw,
                                    amrex::Gpu::HostVector<amrex::Real>& h_avg_vv  , amrex::Gpu::HostVector<amrex::Real>& h_avg_vw,
                                    amrex::Gpu::HostVector<amrex::Real>& h_avg_ww,
                                    amrex::Gpu::HostVector<amrex::Real>& h_avg_uth , amrex::Gpu::HostVector<amrex::Real>& h_avg_vth,
                                    amrex::Gpu::HostVector<amrex::Real>& h_avg_wth, amrex::Gpu::HostVector<amrex::Real>& h_avg_thth,
                                    amrex::Gpu::HostVector<amrex::Real>& h_avg_ku, amrex::Gpu::HostVector<amrex::Real>& h_avg_kv,
                                    amrex::Gpu::HostVector<amrex::Real>& h_avg_kw,
                                    amrex::Gpu::HostVector<amrex::Real>& h_avg_p,
                                    amrex::Gpu::HostVector<amrex::Real>& h_avg_pu, amrex::Gpu::HostVector<amrex::Real>& h_avg_pv,
                                    amrex::Gpu::HostVector<amrex::Real>& h_avg_pw, amrex::Gpu::HostVector<amrex::Real>& h_avg_wthv);

    void derive_stress_profiles (amrex::Gpu::HostVector<amrex::Real>& h_avg_tau11, amrex::Gpu::HostVector<amrex::Real>& h_avg_tau12,
                                 amrex::Gpu::HostVector<amrex::Real>& h_avg_tau13, amrex::Gpu::HostVector<amrex::Real>& h_avg_tau22,
                                 amrex::Gpu::HostVector<amrex::Real>& h_avg_tau23, amrex::Gpu::HostVector<amrex::Real>& h_avg_tau33,
                                 amrex::Gpu::HostVector<amrex::Real>& h_avg_hfx3,  amrex::Gpu::HostVector<amrex::Real>& h_avg_q1fx3,
                                 amrex::Gpu::HostVector<amrex::Real>& h_avg_q2fx3, amrex::Gpu::HostVector<amrex::Real>& h_avg_diss);
    void derive_stress_profiles_stag (amrex::Gpu::HostVector<amrex::Real>& h_avg_tau11, amrex::Gpu::HostVector<amrex::Real>& h_avg_tau12,
                                      amrex::Gpu::HostVector<amrex::Real>& h_avg_tau13, amrex::Gpu::HostVector<amrex::Real>& h_avg_tau22,
                                      amrex::Gpu::HostVector<amrex::Real>& h_avg_tau23, amrex::Gpu::HostVector<amrex::Real>& h_avg_tau33,
                                      amrex::Gpu::HostVector<amrex::Real>& h_avg_hfx3,  amrex::Gpu::HostVector<amrex::Real>& h_avg_q1fx3,
                                      amrex::Gpu::HostVector<amrex::Real>& h_avg_q2fx3, amrex::Gpu::HostVector<amrex::Real>& h_avg_diss);

    // Perform the volume-weighted sum
    amrex::Real
    volWgtSumMF (int lev, const amrex::MultiFab& mf, int comp,
                 const amrex::MultiFab& dJ, const amrex::MultiFab& mfx, const amrex::MultiFab& mfy,
                 bool finemask, bool local = true);

    void
    volWgtColumnSum (int lev, const amrex::MultiFab& mf, int comp,
                     amrex::MultiFab& mf_2d, const amrex::MultiFab& dJ);

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

#ifdef ERF_USE_WW3_COUPLING
    //amrex::Print() <<  " About to call send_to_ww3 from ERF.H" << std::endl;
    void send_to_ww3(int lev);
    //amrex::Print() <<  " About to call read_waves from ERF.H" << std::endl;
    void read_waves(int lev);
    //void send_to_ww3 (int lev);
    //void read_waves (int lev);

    //void send_to_ww3 (int lev);
#endif

    // Interface for advancing the data at one level by one "slow" timestep
    void advance_dycore (int level,
                         amrex::Vector<amrex::MultiFab>& state_old,
                         amrex::Vector<amrex::MultiFab>& state_new,
                         amrex::MultiFab& xvel_old,  amrex::MultiFab& yvel_old,  amrex::MultiFab& zvel_old,
                         amrex::MultiFab& xvel_new,  amrex::MultiFab& yvel_new,  amrex::MultiFab& zvel_new,
                         amrex::MultiFab& source,    amrex::MultiFab& xmom_src,
                         amrex::MultiFab& ymom_src,  amrex::MultiFab& zmom_src,
                         amrex::MultiFab& buoyancy,  amrex::Geometry fine_geom,
                         amrex::Real dt, amrex::Real time);

    void advance_microphysics (int lev,
                               amrex::MultiFab& cons_in,
                               const amrex::Real& dt_advance,
                               const int& iteration,
                               const amrex::Real& time);

    void advance_lsm (int lev,
                      amrex::MultiFab& cons_in,
                      amrex::MultiFab& xvel_in,
                      amrex::MultiFab& yvel_in,
                      const amrex::Real& dt_advance);

    void advance_radiation (int lev,
                            amrex::MultiFab& cons_in,
                            const amrex::Real& dt_advance);

#ifdef ERF_USE_SHOC
    void compute_shoc_tendencies (int lev,
                                  amrex::MultiFab* cons,
                                  amrex::MultiFab* xvel,
                                  amrex::MultiFab* yvel,
                                  amrex::MultiFab* zvel,
                                  amrex::Real* w_subsid,
                                  amrex::MultiFab* tau13,
                                  amrex::MultiFab* tau23,
                                  amrex::MultiFab* hfx3,
                                  amrex::MultiFab* qfx3,
                                  amrex::MultiFab* eddyDiffs,
                                  amrex::MultiFab* z_phys_nd,
                                  const amrex::Real& dt_advance);
#endif

#ifdef ERF_USE_P3
    void compute_p3_tendencies (int lev,
                                amrex::MultiFab& cons_in,
                                amrex::MultiFab& source,
                                const amrex::Real& dt_advance);
#endif

    void build_fine_mask (int lev, amrex::MultiFab& fine_mask);

    void MakeHorizontalAverages ();
    void MakeDiagnosticAverage (amrex::Vector<amrex::Real>& h_havg, amrex::MultiFab& S, int n);
    void derive_upwp (amrex::Vector<amrex::Real>& h_havg);

    // Write plotfile to disk
    void Write3DPlotFile  (int which, PlotFileType plotfile_type, amrex::Vector<std::string> plot_var_names);
    void Write2DPlotFile  (int which, PlotFileType plotfile_type, amrex::Vector<std::string> plot_var_names);

    // Write subvolume of data into a different "plotfile"
    void WriteSubvolume (int isub, amrex::Vector<std::string> subvol_var_names);

    void WriteMultiLevelPlotfileWithTerrain (const std::string &plotfilename,
                                             int nlevels,
                                             const amrex::Vector<const amrex::MultiFab*> &mf,
                                             const amrex::Vector<const amrex::MultiFab*> &mf_nd,
                                             const amrex::Vector<std::string> &varnames,
                                             const amrex::Vector<amrex::Geometry>& my_geom,
                                             amrex::Real time,
                                             const amrex::Vector<int> &level_steps,
                                             const amrex::Vector<amrex::IntVect>& my_ref_ratio,
                                             const std::string &versionName = "HyperCLaw-V1.1",
                                             const std::string &levelPrefix = "Level_",
                                             const std::string &mfPrefix = "Cell",
                                             const amrex::Vector<std::string>& extra_dirs = amrex::Vector<std::string>()) const;


    void WriteGenericPlotfileHeaderWithTerrain (std::ostream &HeaderFile,
                                                int nlevels,
                                                const amrex::Vector<amrex::BoxArray> &bArray,
                                                const amrex::Vector<std::string> &varnames,
                                                const amrex::Vector<amrex::Geometry>& my_geom,
                                                amrex::Real time,
                                                const amrex::Vector<int> &level_steps,
                                                const amrex::Vector<amrex::IntVect>& my_ref_ratio,
                                                const std::string &versionName,
                                                const std::string &levelPrefix,
                                                const std::string &mfPrefix) const;

    void erf_enforce_hse (int lev,
                          amrex::MultiFab& dens, amrex::MultiFab& pres, amrex::MultiFab& pi,
                          amrex::MultiFab& th, amrex::MultiFab& qv,
                          std::unique_ptr<amrex::MultiFab>& z_cc);

#ifdef ERF_USE_NETCDF
    //! Write a timestep to 1D vertical column output for coupling
    void writeToNCColumnFile (int lev,
                              const std::string& colfile_name, amrex::Real xloc, amrex::Real yloc,
                              amrex::Real time);
#endif //ERF_USE_NETCDF

    void init_from_input_sounding (int lev);

    void init_immersed_forcing (int lev);

    void input_sponge (int lev);

    void init_from_hse (int lev);

    void init_thin_body (int lev, const amrex::BoxArray& ba, const amrex::DistributionMapping& dm);

    void FillForecastStateMultiFabs(const int lev,
                                    const std::string& filename,
                                    const std::unique_ptr<amrex::MultiFab>& z_phys_nd,
                                    amrex::Vector<amrex::Vector<amrex::MultiFab>>& weather_forecast_data);

    void WeatherDataInterpolation(const int nlevs,
                                  const amrex::Real time,
                                  amrex::Vector<std::unique_ptr<amrex::MultiFab>>& z_phys_nd,
                                  bool regrid_forces_file_read);

#ifdef ERF_USE_MULTIBLOCK
    // constructor used when ERF is created by a multiblock driver
    // calls AmrCore constructor -> AmrMesh constructor
    ERF (const amrex::RealBox& rb, int max_level_in,
         const amrex::Vector<int>& n_cell_in, int coord,
         const amrex::Vector<amrex::IntVect>& ref_ratio,
         const amrex::Array<int,AMREX_SPACEDIM>& is_per,
         std::string prefix);

    // Advance a block specified number of time steps
    void Evolve_MB (int MBstep, int max_block_step);

    // get the current time values and dt
    amrex::Real get_t_old (int lev) { return t_old[lev]; }
    amrex::Real get_t_new (int lev) { return t_new[lev]; }
    amrex::Real get_dt    (int lev) { return dt[lev];    }

    // Set parmparse prefix for MultiBlock
    void SetParmParsePrefix (std::string name) { pp_prefix = name; }

    // Set 'this' from multiblock container
    void SetMultiBlockPointer (MultiBlockContainer *mbc) { m_mbc = mbc; }

    // Public data copy for MB
    std::vector<amrex::Box> domain_p;

    // Primary solution data data containers
    amrex::Vector<amrex::Vector<amrex::MultiFab> > vars_new;
    amrex::Vector<amrex::Vector<amrex::MultiFab> > vars_old;

    // Velocity time averaged field
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vel_t_avg;
    amrex::Vector<amrex::Real> t_avg_cnt;
#endif

    std::string pp_prefix {"erf"};

    void fill_from_bndryregs (const amrex::Vector<amrex::MultiFab*>& mfs,
                              amrex::Real time);

#ifdef ERF_USE_NETCDF
    // Version in which we always set the velocity values but don't set theta
    void fill_from_realbdy (const amrex::Vector<amrex::MultiFab*>& mfs,
                            amrex::Real time,
                            bool cons_only,
                            int icomp_cons,
                            int ncomp_cons,
                            amrex::IntVect ngvect_cons,
                            amrex::IntVect ngvect_vels);

    // Version in which we only set the boundary value if inflowing, includes theta
    void fill_from_realbdy_upwind (const amrex::Vector<amrex::MultiFab*>& mfs,
                                   amrex::Real time,
                                   bool cons_only,
                                   int icomp_cons,
                                   int ncomp_cons,
                                   amrex::IntVect ngvect_cons,
                                   amrex::IntVect ngvect_vels);
#endif

#ifdef ERF_USE_NETCDF
    void init_from_wrfinput (int lev,
                             amrex::MultiFab& mf_C1H,
                             amrex::MultiFab& mf_C2H,
                             amrex::MultiFab& mf_MUB,
                             amrex::MultiFab& mf_PSFC);
    void init_from_metgrid  (int lev);
    void init_from_ncfile   (int lev);
#endif // ERF_USE_NETCDF

#ifdef ERF_USE_WINDFARM
    void init_windfarm(int lev);
    void advance_windfarm (const amrex::Geometry& a_geom,
                       const amrex::Real& dt_advance,
                       amrex::MultiFab& cons_in,
                       amrex::MultiFab& U_old, amrex::MultiFab& V_old, amrex::MultiFab& W_old,
                       amrex::MultiFab& mf_vars_windfarm,
                       const amrex::MultiFab& mf_Nturb,
                       const amrex::MultiFab& mf_SMark,
                       const amrex::Real& time);
#endif

    void MakeEBGeometry ();
    void make_eb_box ();
    void make_eb_regular ();

    // more flexible version of AverageDown() that lets you average down across multiple levels
    void AverageDownTo (int crse_lev, int scomp, int ncomp); // NOLINT

    // write checkpoint file to disk
    void WriteCheckpointFile () const;

    // read checkpoint file from disk
    void ReadCheckpointFile ();

    // read checkpoint file from disk -- called after instantiating m_SurfaceLayer
    void ReadCheckpointFileSurfaceLayer ();

    void init_zphys            (int lev, amrex::Real time);
    void remake_zphys          (int lev, amrex::Real time, std::unique_ptr<amrex::MultiFab>& temp_zphys_nd);
    void update_terrain_arrays (int lev);

private:

    ///////////////////////////
    // private member functions
    ///////////////////////////

    // read in some parameters from inputs file
    void ReadParameters ();
    void ParameterSanityChecks ();

    // set covered coarse cells/faces to be the average of overlying fine cells/faces
    void AverageDown ();

    void update_diffusive_arrays (int lev, const amrex::BoxArray& ba, const amrex::DistributionMapping& dm);

    void Construct_ERFFillPatchers (int lev);

    void Define_ERFFillPatchers (int lev);

    void init1DArrays ();

    void init_bcs ();
    void init_phys_bcs (bool& rho_read, bool& read_prim_theta);

    void init_custom (int lev);

    void init_uniform (int lev);

    void init_stuff (int lev, const amrex::BoxArray& ba, const amrex::DistributionMapping& dm,
                     amrex::Vector<amrex::MultiFab>& lev_new, amrex::Vector<amrex::MultiFab>& lev_old,
                     amrex::MultiFab& tmp_base_state,
                     std::unique_ptr<amrex::MultiFab>& tmp_zphys_nd);

    // Initialize the Turbulent perturbation
    void turbPert_update (const int lev, const amrex::Real dt);
    void turbPert_amplitude (const int lev);

    // Initialize the integrator object
    void initialize_integrator (int lev, amrex::MultiFab& cons_mf, amrex::MultiFab& vel_mf);

    // Create the physbcs objects
    void make_physbcs (int lev);

    // Initialize the microphysics object
    void initializeMicrophysics (const int&);

#ifdef ERF_USE_WINDFARM
    // Initialize the windfarm object
    void initializeWindFarm (const int&);
#endif

    // Compute a vector of new MultiFabs by copying from valid region and filling ghost cells
    //
    // NOTE: FillPatch takes in an empty MF, and returns cell-centered + velocities (not momenta)
    //
    // This one works only at level = 0 (base state does not change)
    void FillPatchCrseLevel (int lev, amrex::Real time,
                             const amrex::Vector<amrex::MultiFab*>& mfs_vel,
                             bool cons_only=false);

    // This one works only at level > 0 (base state does change)
    void FillPatchFineLevel (int lev, amrex::Real time,
                             const amrex::Vector<amrex::MultiFab*>& mfs_vel,
                             const amrex::Vector<amrex::MultiFab*>& mfs_mom,
                             const amrex::MultiFab& old_base_state,
                             const amrex::MultiFab& new_base_state,
                             bool fillset=true, bool cons_only=false);

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
                                int ng_cons, int ng_vel, bool cons_only, int icomp_cons, int ncomp_cons);

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
    void ComputeDt (int step = -1);

    // get plotfile name
    [[nodiscard]] std::string PlotFileName (int lev) const;

    // set plotfile variable names
    static amrex::Vector<std::string> PlotFileVarNames (amrex::Vector<std::string> plot_var_names);

    // set which variables and derived quantities go into plotfiles
    void setPlotVariables (const std::string& pp_plot_var_names, amrex::Vector<std::string>& plot_var_names);
    void setPlotVariables2D (const std::string& pp_plot_var_names, amrex::Vector<std::string>& plot_var_names);
    // append variables to plot
    void appendPlotVariables (const std::string& pp_plot_var_names, amrex::Vector<std::string>& plot_var_names);

    void setSubVolVariables (const std::string& pp_subvol_var_names, amrex::Vector<std::string>& subvol_var_names);

#ifdef ERF_USE_NETCDF
    //! Create 1D vertical column output for coupling
    void createNCColumnFile (int lev,
                             const std::string& colfile_name, amrex::Real xloc, amrex::Real yloc);

    // Copy from the NC*fabs into the MultiFabs holding the boundary data
    void init_from_wrfbdy (amrex::Vector<amrex::FArrayBox*> x_vel_lateral,
                           amrex::Vector<amrex::FArrayBox*> y_vel_lateral,
                           amrex::Vector<amrex::FArrayBox*> z_vel_lateral,
                           amrex::Vector<amrex::FArrayBox*> T_lateral);

    static amrex::Real start_bdy_time;
    static amrex::Real start_low_time;

    static amrex::Real bdy_time_interval;
    static amrex::Real low_time_interval;

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

    amrex::Vector<amrex::Vector<amrex::FArrayBox>> low_data_zlo;

    // Maximum value of terrain read in from wrfinput at level 0
    amrex::Real z_top = 0.;

#endif // ERF_USE_NETCDF

    amrex::Vector<std::unique_ptr<amrex::MultiFab>> lat_m, lon_m;       // Latitude and Longitude on the grid

    amrex::Vector<std::unique_ptr<amrex::MultiFab>> sinPhi_m, cosPhi_m; // Coriolis factors on the grid

    // Struct for working with the sounding data we take as an input
    InputSoundingData input_sounding_data;

    // Struct for working with the sponge data we take as an input
    InputSpongeData input_sponge_data;

    // Vector (6 planes) of DeviceVectors (ncell in plane) for Dirichlet BC data
    amrex::Vector<amrex::Gpu::DeviceVector<amrex::Real>> xvel_bc_data;
    amrex::Vector<amrex::Gpu::DeviceVector<amrex::Real>> yvel_bc_data;
    amrex::Vector<amrex::Gpu::DeviceVector<amrex::Real>> zvel_bc_data;
    amrex::Vector<amrex::Gpu::DeviceVector<amrex::Real>> th_bc_data;

    // Function to read and populate above vectors (if input file exists)
    void init_Dirichlet_bc_data (const std::string input_file);

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

#ifndef ERF_USE_MULTIBLOCK
    // Array of multifabs to store the solution at each level of refinement;
    // after advancing a level we use "swap".
    amrex::Vector<amrex::Vector<amrex::MultiFab> > vars_new;
    amrex::Vector<amrex::Vector<amrex::MultiFab> > vars_old;

    // Pressure gradient at each level
    amrex::Vector<amrex::Vector<amrex::MultiFab> > gradp;

    // Velocity time averaged field
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vel_t_avg;
    amrex::Vector<amrex::Real> t_avg_cnt;
#endif
    amrex::Vector<std::unique_ptr<MRISplitIntegrator<amrex::Vector<amrex::MultiFab> > > > mri_integrator_mem;

    // Used for anelastic or when we impose an initial projection
    amrex::Vector<amrex::MultiFab> pp_inc;

    // Used only for fast substepping
    amrex::Vector<amrex::MultiFab> lagged_delta_rt;
    amrex::Vector<amrex::MultiFab> avg_xmom;
    amrex::Vector<amrex::MultiFab> avg_ymom;
    amrex::Vector<amrex::MultiFab> avg_zmom;

    // Vector over levels of routines to impose physical boundary conditions
    amrex::Vector<std::unique_ptr<ERFPhysBCFunct_cons>> physbcs_cons;
    amrex::Vector<std::unique_ptr<ERFPhysBCFunct_u>> physbcs_u;
    amrex::Vector<std::unique_ptr<ERFPhysBCFunct_v>> physbcs_v;
    amrex::Vector<std::unique_ptr<ERFPhysBCFunct_w>> physbcs_w;
    amrex::Vector<std::unique_ptr<ERFPhysBCFunct_base>> physbcs_base;

    // Store primitive variable for MOST BC
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> Theta_prim;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> Qv_prim;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> Qr_prim;

    // Scratch space for time integrator
    amrex::Vector<amrex::MultiFab> rU_old;
    amrex::Vector<amrex::MultiFab> rU_new;
    amrex::Vector<amrex::MultiFab> rV_old;
    amrex::Vector<amrex::MultiFab> rV_new;
    amrex::Vector<amrex::MultiFab> rW_old;
    amrex::Vector<amrex::MultiFab> rW_new;

    // amrex::Vector<amrex::MultiFab> xmom_crse_rhs;
    // amrex::Vector<amrex::MultiFab> ymom_crse_rhs;
    amrex::Vector<amrex::MultiFab> zmom_crse_rhs;

    std::unique_ptr<Microphysics> micro;
    amrex::Vector<amrex::Vector<amrex::MultiFab*>> qmoist; // (lev,ncomp) rain_accum, snow_accum, graup_accum

    // Variables for wind farm parametrization models

#ifdef ERF_USE_WINDFARM
    std::unique_ptr<WindFarm> windfarm;
    amrex::Vector<amrex::MultiFab> Nturb;
    amrex::Vector<amrex::MultiFab> vars_windfarm; // Fitch: Vabs, Vabsdt, dudt, dvdt, dTKEdt
                                                  // EWP:                 dudt, dvdt, dTKEdt

    amrex::Vector<amrex::MultiFab> SMark; // A multifab that holds an integer corresponding
                                                   // to the number of the wind turbine to sample
                                                   // velocity upstream of the turbine
#endif

    LandSurface lsm;
    amrex::Vector<std::string> lsm_data_name;                // (ncomp)
    amrex::Vector<amrex::Vector<amrex::MultiFab*>> lsm_data; // (lev,ncomp)
    amrex::Vector<std::string> lsm_flux_name;                // (ncomp)
    amrex::Vector<amrex::Vector<amrex::MultiFab*>> lsm_flux; // (lev,ncomp)

    amrex::Vector<std::unique_ptr<IRadiation>> rad; // Radiation model at each level
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> qheating_rates;  // radiation heating rate source terms (SW, LW)
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> rad_fluxes;      // radiation fluxes (SW up/dn, LW up/dn)
#ifdef ERF_USE_SHOC
    amrex::Vector<std::unique_ptr<SHOCInterface>> shoc_interface; // SHOC model at each level
#endif

#ifdef ERF_USE_P3
    amrex::Vector<std::unique_ptr<P3Interface>> p3_interface; // P3 model at each level
#endif

    // Containers for additional SLM inputs
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> sw_lw_fluxes; // Direct SW (visible, NIR), Diffuse SW (visible, NIR), LW flux
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> solar_zenith; // Solar zenith angle

    bool plot_rad = false;
    int rad_datalog_int = -1;

    // Fillpatcher classes for coarse-fine boundaries
    int cf_width{0};
    int cf_set_width{0};
    amrex::Vector<ERFFillPatcher> FPr_c;
    amrex::Vector<ERFFillPatcher> FPr_u;
    amrex::Vector<ERFFillPatcher> FPr_v;
    amrex::Vector<ERFFillPatcher> FPr_w;

    // Diffusive stresses and Smag
    amrex::Vector<amrex::Vector<std::unique_ptr<amrex::MultiFab>>> Tau;
    amrex::Vector<amrex::Vector<std::unique_ptr<amrex::MultiFab>>> Tau_corr; // correction to tau for implicit
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> eddyDiffs_lev;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> SmnSmn_lev;

    // Sea Surface Temps, Skin Temperature and Land Masks (lev, ntimes)
    amrex::Vector<amrex::Vector<std::unique_ptr<amrex::MultiFab>>>  sst_lev;
    amrex::Vector<amrex::Vector<std::unique_ptr<amrex::MultiFab>>>  tsk_lev;
    amrex::Vector<amrex::Vector<std::unique_ptr<amrex::iMultiFab>>> lmask_lev;

    // Land grid types and urban fraction
    amrex::Vector<amrex::Vector<std::unique_ptr<amrex::iMultiFab>>> land_type_lev;
    amrex::Vector<amrex::Vector<std::unique_ptr<amrex::iMultiFab>>> soil_type_lev;
    amrex::Vector<amrex::Vector<std::unique_ptr<amrex::MultiFab>>>  urb_frac_lev;

    // Other SFS terms
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> SFS_hfx1_lev, SFS_hfx2_lev, SFS_hfx3_lev;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> SFS_diss_lev;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> SFS_q1fx1_lev, SFS_q1fx2_lev, SFS_q1fx3_lev;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> SFS_q2fx3_lev;

    // Grid stretching
    amrex::Vector<amrex::Vector<amrex::Real>> zlevels_stag; // nominal height levels

    // Data structures for terrain-fitted coordinates
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> z_phys_nd;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> z_phys_cc;

    amrex::Vector<std::unique_ptr<amrex::MultiFab>>   detJ_cc;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>>   ax;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>>   ay;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>>   az;

    amrex::Vector<std::unique_ptr<amrex::MultiFab>> z_phys_nd_src;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> z_phys_cc_src;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>>   detJ_cc_src;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>>   ax_src;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>>   ay_src;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>>   az_src;

    amrex::Vector<std::unique_ptr<amrex::MultiFab>> z_phys_nd_new;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>>   detJ_cc_new;

    amrex::Vector<std::unique_ptr<amrex::MultiFab>> z_t_rk;

    // Data structures for terrain immersed forcing
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> terrain_blanking;

    // Wall distance function
    amrex::Vector<std::unique_ptr<amrex::MultiFab>>  walldist;

    // Map scale factors -- vector across levels of vector across types
    amrex::Vector<amrex::Vector<std::unique_ptr<amrex::MultiFab>>> mapfac;

    // Mask to be used to identify cells covered by finer level (if level > 0)
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> fine_mask;

    amrex::Vector<amrex::Vector<amrex::Real>> stretched_dz_h;
    amrex::Vector<amrex::Gpu::DeviceVector<amrex::Real>> stretched_dz_d;

    amrex::Vector<amrex::MultiFab> base_state;
    amrex::Vector<amrex::MultiFab> base_state_new;

    // Wave coupling data
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> Hwave;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> Lwave;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> Hwave_onegrid;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> Lwave_onegrid;
    bool finished_wave = false;

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
    amrex::Array<amrex::Array<amrex::Real, AMREX_SPACEDIM*2>, AMREX_SPACEDIM+NBCVAR_max> m_bc_extdir_vals;

    // These hold the Neumann values at walls which need them ...
    amrex::Array<amrex::Array<amrex::Real, AMREX_SPACEDIM*2>, AMREX_SPACEDIM+NBCVAR_max> m_bc_neumann_vals;

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

    static int last_plot3d_file_step_1;
    static int last_plot3d_file_step_2;
    static int last_plot2d_file_step_1;
    static int last_plot2d_file_step_2;
    static int last_check_file_step;

    static amrex::Real last_plot3d_file_time_1;
    static amrex::Real last_plot3d_file_time_2;
    static amrex::Real last_plot2d_file_time_1;
    static amrex::Real last_plot2d_file_time_2;
    static amrex::Real last_check_file_time;

    amrex::Vector<int> last_subvol_step;
    amrex::Vector<amrex::Real> last_subvol_time;

    static bool plot_file_on_restart;

    // Output control
    const int datwidth = 14;
    const int datprecision = 6;
    const int timeprecision = 13; // e.g., 1-yr LES: 31,536,000 s with dt ~ 0.01 ==> min prec = 10
                                  // epoch time for 2030-01-01 00:00:00 is 1,893,456,000 s with dt ~ 0.01 ==> min prec = 12

    ////////////////
    // runtime parameters

    // Maximum number of steps
    int max_step = -1;

    // Start and stop times
    static amrex::Real start_time;
    static amrex::Real  stop_time;

    bool use_datetime = false;
    const std::string datetime_format = "%Y-%m-%d %H:%M:%S"; // ISO 8601 standard

    // if >= 0 we restart from a checkpoint
    std::string restart_chkfile = "";

    // Time step controls
    static amrex::Real cfl;
    static amrex::Real sub_cfl;
    static amrex::Real init_shrink;
    static amrex::Real change_max;
    static amrex::Real dt_max_initial;
    static amrex::Real dt_max;

    // Fixed dt at each level (only used if positive)
    amrex::Vector<amrex::Real> fixed_dt;
    amrex::Vector<amrex::Real> fixed_fast_dt;
    static int fixed_mri_dt_ratio;

    // how often each level regrids the higher levels of refinement
    // (after a level advances that many time steps)
    int regrid_int = -1;

    // Should we regrid level 0 after restart if Ngrids < Nprocs or
    //    max_grid_size has changed?
    bool regrid_level_0_on_restart = false;

    // plotfile prefix and frequency
    std::string plot3d_file_1 {"plt_1_"};
    std::string plot3d_file_2 {"plt_2_"};
    std::string plot2d_file_1 {"plt2d_1_"};
    std::string plot2d_file_2 {"plt2d_2_"};
    std::string subvol_file {"subvol"};
    bool m_expand_plotvars_to_unif_rr = false;
    int m_plot3d_int_1 = -1;
    int m_plot3d_int_2 = -1;
    int m_plot2d_int_1 = -1;
    int m_plot2d_int_2 = -1;

    amrex::Vector<int> m_subvol_int;
    amrex::Vector<amrex::Real> m_subvol_per;

    amrex::Real m_plot3d_per_1 = -1.0;
    amrex::Real m_plot3d_per_2 = -1.0;
    amrex::Real m_plot2d_per_1 = -1.0;
    amrex::Real m_plot2d_per_2 = -1.0;
    bool m_plot_face_vels = false;

    bool plot_lsm = false;

    // other sampling output control
    int profile_int = -1;
    bool destag_profiles = true;

    // Checkpoint type, prefix and frequency
    std::string check_file {"chk"};
    int m_check_int = -1;
    amrex::Real m_check_per = -1.0;

    amrex::Vector<std::string> subvol3d_var_names;

    amrex::Vector<std::string> plot3d_var_names_1;
    amrex::Vector<std::string> plot3d_var_names_2;
    amrex::Vector<std::string> plot2d_var_names_1;
    amrex::Vector<std::string> plot2d_var_names_2;
    const amrex::Vector<std::string> cons_names     {"density", "rhotheta", "rhoKE", "rhoadv_0",
                                                     "rhoQ1", "rhoQ2", "rhoQ3",
                                                     "rhoQ4", "rhoQ5", "rhoQ6"};

    // **************************************************************************************
    // NOTE: The order of variable names here **MUST MATCH THE ORDER** in IO/ERF_Plotfile.cpp
    // **************************************************************************************
    const amrex::Vector<std::string> derived_names {"soundspeed", "temp", "theta", "KE", "scalar",
                                                    "vorticity_x","vorticity_y","vorticity_z",
                                                    "magvel", "divU",
                                                    "pres_hse", "dens_hse", "theta_hse", "pressure", "pert_pres", "pert_dens",
                                                    "eq_pot_temp",
#ifdef ERF_USE_WINDFARM
                                                    "num_turb", "SMark0", "SMark1",
#endif
                                                    "dpdx", "dpdy", "dpdz", "pres_hse_x", "pres_hse_y",
                                                    "z_phys", "detJ", "mapfac", "lat_m", "lon_m",
                                                    // Time averaged velocity
                                                    "u_t_avg", "v_t_avg", "w_t_avg", "umag_t_avg",
                                                    // eddy viscosity
                                                    "nut",
                                                    // eddy diffusivity of momentum
                                                    "Kmv","Kmh",
                                                    // eddy diffusivity of heat
                                                    "Khv","Khh",
                                                    // turbulence lengthscale
                                                    "Lturb",
                                                    // wall distance
                                                    "walldist",
                                                    // dissipation
                                                    "diss",
                                                    // moisture vars
                                                    "moist_density", "qv", "qc", "qi", "qrain", "qsnow", "qgraup", "qt", "qn", "qp", "qsat",
                                                    "rain_accum", "snow_accum", "graup_accum", "reflectivity", "max_reflectivity",
                                                    // Terrain IB mask
                                                    "terrain_IB_mask",
                                                    // EB variables
                                                    "volfrac"
#ifdef ERF_COMPUTE_ERROR
                                                    // error vars
                                                    ,"xvel_err", "yvel_err", "zvel_err", "pp_err"
#endif
                                                    ,"qsrc_sw", "qsrc_lw"
                                                   };

    // **************************************************************************************
    // NOTE: The order of variable names here **MUST MATCH THE ORDER** in IO/ERF_Plotfile.cpp
    // **************************************************************************************
    const amrex::Vector<std::string> derived_names_2d {
        "z_surf", "landmask", "mapfac", "lat_m", "lon_m",
        "u_star", "w_star", "t_star", "q_star", "Olen", "pblh",
        "t_surf", "q_surf", "z0", "OLR", "sens_flux", "laten_flux",
        "surf_pres", "integrated_qv"
    };

    // **************************************************************************************
    // NOTE: The order of variable names here **MUST MATCH THE ORDER** in IO/ERF_WriteSubVolume.cpp
    // **************************************************************************************
    const amrex::Vector<std::string> derived_subvol_names {"soundspeed", "temp", "theta", "KE", "scalar"};

    // algorithm choices
    static SolverChoice solverChoice;

    // Turbulent perturbation structure
    TurbulentPerturbation turbPert;

#ifdef ERF_USE_PARTICLES
    // Particle container with all particle species
    ParticleData particleData;

    // variables and functions for tracers particles
    bool  m_use_tracer_particles; /*!< tracer particles that advect with flow with optional sedimentation */

    /*! Read tracer and hydro particles parameters */
    void readTracersParams();

    /*! Initialize tracer and hydro particles */
    void initializeTracers ( amrex::ParGDBBase*,
                             const amrex::Vector<std::unique_ptr<amrex::MultiFab>>&,
                             const amrex::Real time);

    /*! Restart tracer and hydro particles */
    void restartTracers ( amrex::ParGDBBase*, const std::string& );

    /*! Evolve tracers and hydro particles */
    void evolveTracers( int,
                        amrex::Real,
                        amrex::Vector<amrex::Vector<amrex::MultiFab>>&,
                        const amrex::Vector<std::unique_ptr<amrex::MultiFab>>& );

#endif

#ifdef ERF_USE_MULTIBLOCK
    MultiBlockContainer *m_mbc = nullptr;
#endif

    static int verbose;
    static int mg_verbose;
    static bool use_fft;

    static int check_for_nans;

    // Diagnostic output interval
    static int sum_interval;
    static int pert_interval;
    static amrex::Real sum_per;

    // Write in native AMReX or NetCDF format for each plotfile
    static PlotFileType plotfile3d_type_1;
    static PlotFileType plotfile3d_type_2;
    static PlotFileType plotfile2d_type_1;
    static PlotFileType plotfile2d_type_2;

    int file_name_digits = 5;
    bool use_real_time_in_pltname = false;

    static StateInterpType interpolation_type;

    // sponge_type:  "input_sponge"
    static std::string sponge_type;

    // NetCDF initialization (wrfinput) file
    static amrex::Vector<amrex::Vector<std::string>> nc_init_file;
    static amrex::Vector<amrex::Vector<int>> have_read_nc_init_file;

    // NetCDF initialization (wrfbdy/met_em) file
    static std::string nc_bdy_file;
    int real_width{0};
    int real_set_width{0};
    bool real_extrap_w{true};

    // NetCDF initialization of zlo boundary values
    static std::string nc_low_file;

    // Options for vertical interpolation of met_em*.nc data.
    bool        metgrid_debug_quiescent{false};
    bool        metgrid_debug_isothermal{false};
    bool        metgrid_debug_dry{false};
    bool        metgrid_debug_psfc{false};
    bool        metgrid_debug_msf{false};
    bool        metgrid_interp_theta{false};
    bool        metgrid_basic_linear{false};
    bool        metgrid_use_below_sfc{true};
    bool        metgrid_use_sfc{true};
    bool        metgrid_retain_sfc{false};
    amrex::Real metgrid_proximity{500.0};
    int         metgrid_order{2};
    int         metgrid_force_sfc_k{6};

    amrex::Vector<amrex::BoxArray> ba1d;
    amrex::Vector<amrex::BoxArray> ba2d;
    std::unique_ptr<amrex::MultiFab> mf_C1H;
    std::unique_ptr<amrex::MultiFab> mf_C2H;
    std::unique_ptr<amrex::MultiFab> mf_MUB;

    amrex::Vector<std::unique_ptr<amrex::MultiFab>> mf_PSFC;

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

    // Function to read and populate above host vectors (if input file exists)
    void init_geo_wind_profile (const std::string input_file,
                                amrex::Vector<amrex::Real>& u_geos,
                                amrex::Gpu::DeviceVector<amrex::Real>& u_geos_d,
                                amrex::Vector<amrex::Real>& v_geos,
                                amrex::Gpu::DeviceVector<amrex::Real>& v_geos_d,
                                const amrex::Geometry& lgeom,
                                const amrex::Vector<amrex::Real>& zlev_stag);

    // This is a vector over levels of vectors across quantities of Vectors
    amrex::Vector<amrex::Vector<amrex::Vector<amrex::Real> > > h_rayleigh_ptrs;
    amrex::Vector<amrex::Vector<amrex::Vector<amrex::Real> > > h_sponge_ptrs;

    // These are vectors over levels of Vectors
    amrex::Vector<amrex::Vector<amrex::Real> > h_sinesq_ptrs;
    amrex::Vector<amrex::Vector<amrex::Real> > h_sinesq_stag_ptrs;

    // This is a vector over levels of vectors across quantities of DeviceVectors
    amrex::Vector<amrex::Vector<amrex::Gpu::DeviceVector<amrex::Real> > > d_rayleigh_ptrs;
    amrex::Vector<amrex::Vector<amrex::Gpu::DeviceVector<amrex::Real> > > d_sponge_ptrs;

    // These are vectors over levels of DeviceVectors
    amrex::Vector<amrex::Gpu::DeviceVector<amrex::Real> > d_sinesq_ptrs;
    amrex::Vector<amrex::Gpu::DeviceVector<amrex::Real> > d_sinesq_stag_ptrs;

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
    std::unique_ptr<SurfaceLayer>     m_SurfaceLayer = nullptr;
    amrex::Vector<std::unique_ptr<ForestDrag>> m_forest_drag;

    //
    // Holds info for dynamically generated tagging criteria
    //
    static amrex::Vector<amrex::AMRErrorTag> ref_tags;

    amrex::Vector<amrex::Vector<amrex::BoxArray>> subdomains;

    amrex::Vector<amrex::Real> dz_min;

    static AMREX_FORCE_INLINE
    int
    ComputeGhostCells (const SolverChoice& sc)
    {
        int ngrow = 0;

        if (sc.use_num_diff)
        {
            ngrow = 3;
        } else {
            if (
                 (sc.advChoice.dycore_horiz_adv_type    == AdvType::Centered_6th)
              || (sc.advChoice.dycore_vert_adv_type     == AdvType::Centered_6th)
              || (sc.advChoice.moistscal_horiz_adv_type == AdvType::Centered_6th)
              || (sc.advChoice.moistscal_horiz_adv_type == AdvType::Centered_6th)
              || (sc.advChoice.dryscal_horiz_adv_type   == AdvType::Centered_6th)
              || (sc.advChoice.dryscal_vert_adv_type    == AdvType::Centered_6th) )
            { ngrow = 3; }
            else if (
                 (sc.advChoice.dycore_horiz_adv_type    == AdvType::Upwind_5th)
              || (sc.advChoice.dycore_vert_adv_type     == AdvType::Upwind_5th)
              || (sc.advChoice.moistscal_horiz_adv_type == AdvType::Upwind_5th)
              || (sc.advChoice.moistscal_horiz_adv_type == AdvType::Upwind_5th)
              || (sc.advChoice.dryscal_horiz_adv_type   == AdvType::Upwind_5th)
              || (sc.advChoice.dryscal_vert_adv_type    == AdvType::Upwind_5th) )
            { ngrow = 3; }
            else if (
                  (sc.advChoice.dryscal_horiz_adv_type   == AdvType::Weno_5)
               || (sc.advChoice.dryscal_vert_adv_type    == AdvType::Weno_5)
               || (sc.advChoice.moistscal_horiz_adv_type == AdvType::Weno_5)
               || (sc.advChoice.moistscal_vert_adv_type  == AdvType::Weno_5)
               || (sc.advChoice.moistscal_horiz_adv_type == AdvType::Weno_5Z)
               || (sc.advChoice.moistscal_vert_adv_type  == AdvType::Weno_5Z)
               || (sc.advChoice.dryscal_horiz_adv_type   == AdvType::Weno_5Z)
               || (sc.advChoice.dryscal_vert_adv_type    == AdvType::Weno_5Z) )
            { ngrow = 3; }
            else if (
                  (sc.advChoice.dryscal_horiz_adv_type   == AdvType::Weno_7)
               || (sc.advChoice.dryscal_vert_adv_type    == AdvType::Weno_7)
               || (sc.advChoice.moistscal_horiz_adv_type == AdvType::Weno_7)
               || (sc.advChoice.moistscal_vert_adv_type  == AdvType::Weno_7)
               || (sc.advChoice.moistscal_horiz_adv_type == AdvType::Weno_7Z)
               || (sc.advChoice.moistscal_vert_adv_type  == AdvType::Weno_7Z)
               || (sc.advChoice.dryscal_horiz_adv_type   == AdvType::Weno_7Z)
               || (sc.advChoice.dryscal_vert_adv_type    == AdvType::Weno_7Z) )
            { ngrow = 4; }
            else
            {
                if (sc.terrain_type == TerrainType::EB){
                    ngrow = 3;
                } else {
                    ngrow = 2;
                }
            }
        }

        return ngrow;
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
    std::ostream&
    DerDataLog (int i)
    {
        return *der_datalog[i];
    }

    AMREX_FORCE_INLINE
    int
    NumDataLogs () noexcept
    {
        return datalog.size();
    }

    AMREX_FORCE_INLINE
    int
    NumDerDataLogs () noexcept
    {
        return der_datalog.size();
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

    void setRecordDerDataInfo (int i, const std::string& filename) // NOLINT
    {
        if (amrex::ParallelDescriptor::IOProcessor())
        {
            der_datalog[i] = std::make_unique<std::fstream>();
            der_datalog[i]->open(filename.c_str(),std::ios::out|std::ios::app);
            if (!der_datalog[i]->good()) {
                amrex::FileOpenFailed(filename);
            }
        }
        amrex::ParallelDescriptor::Barrier("ERF::setRecordDerDataInfo");
    }

    void setRecordEnergyDataInfo (int i, const std::string& filename) // NOLINT
    {
        if (amrex::ParallelDescriptor::IOProcessor())
        {
            tot_e_datalog[i] = std::make_unique<std::fstream>();
            tot_e_datalog[i]->open(filename.c_str(),std::ios::out|std::ios::app);
            if (!tot_e_datalog[i]->good()) {
                amrex::FileOpenFailed(filename);
            }
        }
        amrex::ParallelDescriptor::Barrier("ERF::setRecordEnergyDataInfo");
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

    // Data sampler for line and plane output
    int line_sampling_interval = -1;
    int plane_sampling_interval = -1;
    amrex::Real line_sampling_per = -1.0;
    amrex::Real plane_sampling_per = -1.0;
    std::unique_ptr<LineSampler> line_sampler = nullptr;
    std::unique_ptr<PlaneSampler> plane_sampler = nullptr;

    amrex::Vector<std::unique_ptr<std::fstream> >       datalog;
    amrex::Vector<std::unique_ptr<std::fstream> >   der_datalog;
    amrex::Vector<std::unique_ptr<std::fstream> > tot_e_datalog;
    amrex::Vector<std::string>       datalogname;
    amrex::Vector<std::string>   der_datalogname;
    amrex::Vector<std::string> tot_e_datalogname;

    amrex::Vector<std::unique_ptr<std::fstream> > sampleptlog;
    amrex::Vector<std::string> sampleptlogname;
    amrex::Vector<amrex::IntVect> samplepoint;

    amrex::Vector<std::unique_ptr<std::fstream> > samplelinelog;
    amrex::Vector<std::string> samplelinelogname;
    amrex::Vector<amrex::IntVect> sampleline;

    //! The filename of the ith datalog file.
    [[nodiscard]] std::string DataLogName (int i) const noexcept { return datalogname[i]; }
    [[nodiscard]] std::string DerDataLogName (int i) const noexcept { return der_datalogname[i]; }

    //! The filename of the ith sampleptlog file.
    [[nodiscard]] std::string SamplePointLogName (int i) const noexcept { return sampleptlogname[i]; }

    //! The filename of the ith samplelinelog file.
    [[nodiscard]] std::string SampleLineLogName (int i) const noexcept { return samplelinelogname[i]; }

    // array of EB objects
    amrex::Vector<std::unique_ptr<eb_>> eb;

    [[nodiscard]] eb_ const& get_eb (int lev) const noexcept {
        AMREX_ASSERT(lev >= 0 && lev < eb.size() && eb[lev] != nullptr);
        return *eb[lev];
    }

    [[nodiscard]] amrex::EBFArrayBoxFactory const&
    EBFactory (int lev) const noexcept {
        return *(eb[lev]->get_const_factory());
    }

    [[nodiscard]] static int nghost_eb_basic ()
    { return 5; }

    // We need 5 for doing StateRedistribution; otherwise 4 would be enough
    [[nodiscard]] static int nghost_eb_volume ()
    { return 5; }

    [[nodiscard]] static int nghost_eb_full ()
    { return 4; }

#ifdef ERF_USE_FFT
    amrex::Vector<std::unique_ptr<amrex::FFT::Poisson<amrex::MultiFab>>> m_3D_poisson;
    amrex::Vector<std::unique_ptr<amrex::FFT::PoissonHybrid<amrex::MultiFab>>> m_2D_poisson;
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
