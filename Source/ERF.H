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
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_FluxRegister.H>
#include <AMReX_ErrorList.H>

#ifdef AMREX_MEM_PROFILING
#include <AMReX_MemProfiler.H>
#endif

#include <IndexDefines.H>
#include <TimeInterpolatedData.H>
#include <ERF_PhysBCFunct.H>
#include <AMReX_TimeIntegrator.H>
#include <TimeIntegration.H>
#include <DataStruct.H>
#include <ABLFieldInit.H>
#include <Derive.H>
//#include <prob.H>

#include <iostream>

#ifdef AMREX_LAZY
#include <AMReX_Lazy.H>
#endif

using namespace amrex;

namespace InterpType {
    enum {
        PCInterp = 0,
        NodeBilinear,
        CellConservativeLinear,
        CellBilinear,
        CellQuadratic,
        CellConservativeProtected,
        CellConservativeQuartic
    };
}

template<typename V, typename T>
bool containerHasElement(const V& iterable, const T& query) {
    return std::find(iterable.begin(), iterable.end(), query) != iterable.end();
}

class ERF
    : public amrex::AmrCore
{
public:

    ////////////////
    // public member functions

    // constructor - reads in parameters from inputs file
    //             - sizes multilevel arrays and data structures
    ERF ();
    virtual ~ERF();

    // Advance solution to final time
    void Evolve ();

    // Tag cells for refinement
    virtual void ErrorEst (int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow) override;

    // Initialize multilevel data
    void InitData ();

    // Called after every level 0 timestep
    void post_timestep (int nstep, Real time, Real dt_lev);

    // Diagnostics
    void sum_integrated_quantities(Real time);

    // Perform the volume-weighted sum
    amrex::Real
    volWgtSumMF(int lev,
      const amrex::MultiFab& mf, int comp, bool local, bool finemask);

    // Decide if it is time to take an action
    bool is_it_time_for_action(int nstep, amrex::Real time, amrex::Real dt,
                               int action_interval, amrex::Real action_per);

    // Make a new level using provided BoxArray and DistributionMapping and
    // fill with interpolated coarse level data.
    // overrides the pure virtual function in AmrCore
    virtual void MakeNewLevelFromCoarse (int lev, amrex::Real time, const amrex::BoxArray& ba,
                     const amrex::DistributionMapping& dm) override;

    // Remake an existing level using provided BoxArray and DistributionMapping and
    // fill with existing fine and coarse data.
    // overrides the pure virtual function in AmrCore
    virtual void RemakeLevel (int lev, amrex::Real time, const amrex::BoxArray& ba,
                  const amrex::DistributionMapping& dm) override;

    // Delete level data
    // overrides the pure virtual function in AmrCore
    virtual void ClearLevel (int lev) override;

    // Make a new level from scratch using provided BoxArray and DistributionMapping.
    // Only used during initialization.
    // overrides the pure virtual function in AmrCore
    virtual void MakeNewLevelFromScratch (int lev, amrex::Real time, const amrex::BoxArray& ba,
                      const amrex::DistributionMapping& dm) override;

private:

    ////////////////
    // private member functions

    // read in some parameters from inputs file
    void ReadParameters();

    // set covered coarse cells to be the average of overlying fine cells
    void AverageDown ();

    void init1DArrays();

    void init_bcs();

    // more flexible version of AverageDown() that lets you average down across multiple levels
    void AverageDownTo (int crse_lev);

    // compute a new multifab by coping in phi from valid region and filling ghost cells
    // works for single level and 2-level cases (fill fine grid ghost by interpolating from coarse)
    void FillPatch (int lev, amrex::Real time, amrex::MultiFab& mf, int icomp, int ncomp, int var_idx);

    // compute a new multifab by coping in data from valid region and filling ghost cells
    // works for single level and 2-level cases (fill fine grid ghost by interpolating from coarse)
    // unlike FillPatch, FillIntermediatePatch will use the supplied multifab instead of fine level data.
    // This is to support filling boundary cells at an intermediate time between old/new times
    // on the fine level when valid data at a specific time is already available (such as
    // at each RK stage when integrating between initial and final times at a given level).
    void FillIntermediatePatch (int lev, Real time, MultiFab& mf, int icomp, int ncomp, int var_idx);

    // fill an entire multifab by interpolating from the coarser level
    // this comes into play when a new level of refinement appears
    void FillCoarsePatch (int lev, amrex::Real time, amrex::MultiFab& mf, int icomp, int ncomp, int var_idx);

    // fill all multifabs (and all components) in a vector of multifabs corresponding to the
    // grid variables defined in vars_old and vars_new just as FillCoarsePatch.
    void FillCoarsePatchAllVars (int lev, Real time, Vector<MultiFab>& vmf);

    // utility to copy in data from old and/or new state into another multifab
    TimeInterpolatedData GetDataAtTime (int lev, Real time);

    // get a vector of BCRec that we can pass to FillPatch depending on what type of variable we want to fill
    Vector<BCRec> GetBCVarVector (int var_idx);

    // advance a level by dt
    // includes a recursive call for finer levels
    void timeStep (int lev, amrex::Real time, int iteration);

    // advance a single level for a single time step
    void Advance (int lev, amrex::Real time, amrex::Real dt_lev, int iteration, int ncycle);

    void erf_advance(int level,
                      MultiFab& cons_old,  MultiFab& cons_new,
                      MultiFab& xvel_old, MultiFab& yvel_old, MultiFab& zvel_old,
                      MultiFab& xvel_new, MultiFab& yvel_new, MultiFab& zvel_new,
                      MultiFab& xmom_crse, MultiFab& ymom_crse, MultiFab& zmom_crse,
                      MultiFab& source,
                      std::array< MultiFab, AMREX_SPACEDIM>& flux,
                      const amrex::Geometry crse_geom,
                      const amrex::Geometry fine_geom,
                      const amrex::IntVect ref_ratio,
                      const amrex::Real dt, const amrex::Real time,
                      amrex::InterpFaceRegister* ifr,
                      const SolverChoice& solverChoice,
                      const amrex::Real* dptr_dens_hse,
                      const amrex::Real* dptr_pres_hse,
                      const amrex::Real* dptr_rayleigh_tau,
                      const amrex::Real* dptr_rayleigh_ubar,
                      const amrex::Real* dptr_rayleigh_vbar,
                      const amrex::Real* dptr_rayleigh_thetabar);

    //! Initialize HSE
    void initHSE();
    void erf_enforce_hse(int lev,
                         amrex::Vector<amrex::Real>& dens,
                         amrex::Vector<amrex::Real>& pres);

    //! Initialize Rayleigh damping profiles
    void initRayleigh();

    // a wrapper for estTimeStep()
    void ComputeDt ();

    // compute dt from CFL considerations
    Real estTimeStep (int lev) const;

    // get plotfile name
    std::string PlotFileName (int lev) const;

    // set plotfile variables names
    amrex::Vector<std::string> PlotFileVarNames () const;

    // set which variables and derived quantities go into plotfiles
    void setPlotVariables ();

    // write plotfile to disk
    void WritePlotFile () const;

#ifdef ERF_USE_NETCDF
    //! Write plotfile using NETCDF
    void writeNCPlotFile(const std::string& dir,
                         const Vector<const MultiFab*> &mf,
                         const Vector<std::string> &varnames,
                         const Vector<int> level_steps, const Real time) const;

    //! Create 1D vertical column output for coupling
    void createNCColumnFile(const int lev,
                            const std::string& colfile_name, const amrex::Real xloc, const amrex::Real yloc);

    //! Write a timestep to 1D vertical column output for coupling
    void writeToNCColumnFile(const int lev,
                             const std::string& colfile_name, const amrex::Real xloc, const amrex::Real yloc,
                             const amrex::Real time);
#endif

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
    void post_update (amrex::MultiFab& state_mf, const amrex::Real time, const amrex::Geometry& geom);
    void fill_rhs (amrex::MultiFab& rhs_mf, const amrex::MultiFab& state_mf, const amrex::Real time, const amrex::Geometry& geom);

    ////////////////
    // private data members

    amrex::Vector<int> istep;      // which step?
    amrex::Vector<int> nsubsteps;  // how many substeps on each level?

    // keep track of old time, new time, and time step at each level
    amrex::Vector<amrex::Real> t_new;
    amrex::Vector<amrex::Real> t_old;
    amrex::Vector<amrex::Real> dt;

    // array of multifabs to store the solution at each level of refinement
    // after advancing a level we use "swap".
    amrex::Vector<amrex::Vector<amrex::MultiFab> > vars_new;
    amrex::Vector<amrex::Vector<amrex::MultiFab> > vars_old;

    // array of time integrator objects at each level of refinement
    amrex::Vector<std::unique_ptr<TimeIntegrator<Vector<MultiFab> > > > integrator;

    // array of flux registers
    amrex::Vector<amrex::FluxRegister*> flux_registers;

    // A BCRec is essentially a 2*DIM integer array storing the boundary
    // condition type at each lo/hi walls in each direction
    // Here bcs[0]/bcs_d[0] correspond to bcs for velocity
    //      bcs[1]/bcs_d[1] correspond to bcs for cell-centered quantities
    amrex::Vector           <amrex::Vector           <amrex::BCRec>> bcs;
    amrex::Vector<amrex::Gpu::DeviceVector<amrex::BCRec>> bcs_d;

    // These hold the Dirichlet values at walls which need them ...
    amrex::GpuArray<amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>
                                               , AMREX_SPACEDIM*2> m_bc_vels;
    amrex::GpuArray<amrex::GpuArray<amrex::Real, NVAR>
                                               , AMREX_SPACEDIM*2> m_bc_cons;

    enum struct BC {
        symmetry, inflow, outflow, no_slip_wall, slip_wall, periodic, MOST, undefined
    };

    // These are the "physical" boundary condition types (e.g. "inflow")
    amrex::GpuArray<BC, AMREX_SPACEDIM*2> m_bc_type;

    // Interpolater to use between levels
    amrex::Interpolater* mapper = nullptr;

    int last_plot_file_step;
    int last_check_file_step;
    int plot_file_on_restart = 1;

    ////////////////
    // runtime parameters

    // maximum number of steps and stop time
    int max_step = std::numeric_limits<int>::max();
    amrex::Real stop_time = std::numeric_limits<amrex::Real>::max();

    // if >= 0 we restart from a checkpoint
    std::string restart_chkfile = "";

    // Time step controls
    static amrex::Real cfl;
    static amrex::Real init_shrink;
    static amrex::Real change_max;

    // Fixed dt for level 0 timesteps (only used if positive)
    static amrex::Real fixed_dt;

    // how often each level regrids the higher levels of refinement
    // (after a level advances that many time steps)
    int regrid_int = 2;

    // plotfile prefix and frequency
    std::string plot_file {"plt"};
    int plot_int = -1;

    // checkpoint prefix and frequency
    std::string check_file {"chk"};
    int check_int = -1;

    amrex::Vector<std::string> plot_state_names;
    amrex::Vector<std::string> plot_deriv_names;
    const amrex::Vector<std::string> velocity_names {"x_velocity", "y_velocity", "z_velocity"};
    const amrex::Vector<std::string> cons_names {"density", "rhotheta", "rhoKE", "rhoadv_0"};
    const amrex::Vector<std::string> derived_names {"pressure", "soundspeed", "temp", "theta",
                                                    "pres_hse", "dens_hse", "pert_pres", "pert_dens"};

    // interpolation type between levels
    int interpolation_type = 0;

    // algorithm choices
    static SolverChoice solverChoice;

    static int verbose;

    // mesh refinement
    static std::string coupling_type;
    static int do_reflux;
    static int do_avg_down;

    // Diagnostic output interval
    static int sum_interval;
    static amrex::Real sum_per;

    // Native or NetCDF
    static std::string plotfile_type;
    static std::string checkpoint_type;

    // 1D CDF output (for ingestion in AMR-Wind)
    static int         output_1d_column;
    static int         column_interval;
    static amrex::Real column_per;
    static amrex::Real column_loc_x;
    static amrex::Real column_loc_y;
    static std::string column_file_name;

    // Base state quantities: density and pressure, 1 array per AMR level
    static amrex::Vector<amrex::Vector<amrex::Real> > h_dens_hse;
    static amrex::Vector<amrex::Vector<amrex::Real> > h_pres_hse;
    static amrex::Vector<amrex::Gpu::DeviceVector<amrex::Real> > d_dens_hse;
    static amrex::Vector<amrex::Gpu::DeviceVector<amrex::Real> > d_pres_hse;

    static int ng_dens_hse;
    static int ng_pres_hse;

    // Mean quantities and tau for Rayleigh damping
    static amrex::Vector<amrex::Vector<amrex::Real> > h_rayleigh_tau;
    static amrex::Vector<amrex::Vector<amrex::Real> > h_rayleigh_ubar;
    static amrex::Vector<amrex::Vector<amrex::Real> > h_rayleigh_vbar;
    static amrex::Vector<amrex::Vector<amrex::Real> > h_rayleigh_thetabar;
    static amrex::Vector<amrex::Gpu::DeviceVector<amrex::Real> > d_rayleigh_tau;
    static amrex::Vector<amrex::Gpu::DeviceVector<amrex::Real> > d_rayleigh_ubar;
    static amrex::Vector<amrex::Gpu::DeviceVector<amrex::Real> > d_rayleigh_vbar;
    static amrex::Vector<amrex::Gpu::DeviceVector<amrex::Real> > d_rayleigh_thetabar;

    static void refinement_criteria_setup();

    //
    // Holds info for dynamically generated tagging criteria
    //
    static amrex::Vector<amrex::AMRErrorTag> ref_tags;

    //
    // Build a mask that zeroes out values on a coarse level underlying
    //     grids on the next finest level
    //
    amrex::MultiFab& build_fine_mask(int lev);
    amrex::MultiFab fine_mask;

    AMREX_FORCE_INLINE
    int
    ComputeGhostCells(const int& spatial_order) {
      int nGhostCells;

      switch (spatial_order) {
        case 2:
          nGhostCells = 2; // We need this many to compute the eddy viscosity in the ghost cells
          break;
        case 3:
          nGhostCells = 2;
          break;
        case 4:
          nGhostCells = 2;
          break;
        case 5:
          nGhostCells = 2;
          break;
        case 6:
          nGhostCells = 3;
          break;
        default:
          amrex::Error("Must specify spatial order to be 2,3,4,5 or 6");
      }

      return nGhostCells;
    }

    AMREX_FORCE_INLINE
    amrex::FluxRegister&
    get_flux_reg(int lev)
    {
      return *flux_registers[lev];
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

    static amrex::Real startCPUTime;
    static amrex::Real previousCPUTimeUsed;

    amrex::Real
    getCPUTime() const
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

    void setRecordDataInfo (int i, const std::string& filename)
    {
        if (ParallelDescriptor::IOProcessor())
        {
            datalog[i] = std::make_unique<std::fstream>();
            datalog[i]->open(filename.c_str(),std::ios::out|std::ios::app);
            if (!datalog[i]->good()) {
                amrex::FileOpenFailed(filename);
            }
        }
        ParallelDescriptor::Barrier("ERF::setRecordDataInfo");
    }

    Vector<std::unique_ptr<std::fstream> > datalog;
    Vector<std::string> datalogname;

    //! The filename of the ith datalog file.
    const std::string DataLogName (int i) const noexcept { return datalogname[i]; }

public:
    // ABL Initialization (from amr-wind)
    static ABLFieldInit ablinit;

    static bool lo_z_is_dirichlet;
    static bool hi_z_is_dirichlet;
    static bool init_abl;

    void writeJobInfo(const std::string& dir) const;
    static void writeBuildInfo(std::ostream& os);
};

#endif
