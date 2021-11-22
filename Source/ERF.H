/**
 * \file ERF.H
 */

#ifndef _ERF_H_
#define _ERF_H_

#include <iostream>

#include <AMReX_BC_TYPES.H>
#include <AMReX_AmrLevel.H>
#include <AMReX_iMultiFab.H>

#include <AMReX_YAFluxRegister.H>

#ifdef AMREX_LAZY
#include <AMReX_Lazy.H>
#endif

#include "SolverStruct.H"
#include "IndexDefines.H"
#include "IOManager.H"

#include "BC/BCBase.H"
#include "BC/BCDummy.H"
#include "BC/BCInterior.H"
#include "BC/BCDirichlet.H"
#include "BC/BCSlipWall.H"
#include "BC/BCNoSlipWall.H"
#include "BC/BCSimSlipWall.H"

using std::istream;
using std::ostream;

#define UserBC 6

#define ENUM_BUILDER(name, ...) enum name { __VA_ARGS__ }; \
inline std::string EnumToString(name value) {              \
std::string str = #__VA_ARGS__;                            \
std::string rs;                                            \
int len = str.length();                                    \
std::vector<std::string> strings;                          \
std::ostringstream ss;                                     \
for(int i = 0; i < len; i ++) {                            \
    if(isspace(str[i])) continue;                          \
    else if(str[i] == ',') {                               \
        strings.push_back(ss.str());                       \
        ss.str(std::string());                             \
        }                                                  \
    else ss<<str[i];                                       \
}                                                          \
strings.push_back(ss.str());                               \
rs = strings[static_cast<int>(value)];                     \
return rs; }


ENUM_BUILDER(StateType,
             State_Type,
             X_Vel_Type,
             Y_Vel_Type,
             Z_Vel_Type,
             Work_Estimate_Type)

/**
 * AmrLevel-derived class for solving equations for atmospheric flows
 */
class ERF : public amrex::AmrLevel
{
public:
  //! Default constructor.  Builds invalid object.
  ERF();

  static SolverChoice solverChoice;

  //! The basic constructor.
  ERF(
    amrex::Amr& papa,
    int lev,
    const amrex::Geometry& level_geom,
    const amrex::BoxArray& bl,
    const amrex::DistributionMapping& dm,
    amrex::Real time);

  //! The destructor.
  virtual ~ERF();

  //! Restart from a checkpoint file.
  virtual void
  restart(amrex::Amr& papa, istream& is, bool bReadSpecial = false) override;

  //! This is called only when we restart from an old checkpoint.
  virtual void
  set_state_in_checkpoint(amrex::Vector<int>& state_in_checkpoint) override;

  //! Call AmrLevel::checkPoint and then add radiation info
  virtual void checkPoint(
    const std::string& dir,
    std::ostream& os,
    amrex::VisMF::How how,
    bool dump_old) override;

  virtual void setPlotVariables() override;

  //! Write a plotfile to specified directory.
  virtual void writePlotFile(
    const std::string& dir, ostream& os, amrex::VisMF::How how) override;
  virtual void writeSmallPlotFile(
    const std::string& dir, ostream& os, amrex::VisMF::How how) override;
  void writeJobInfo(const std::string& dir);

  IOManager& get_io_manager() { return *io_mgr; }
  const IOManager& get_io_manager() const { return *io_mgr; }

  //! Define data descriptors.
  static void variableSetUp();

  //! Cleanup data descriptors at end of run.
  static void variableCleanUp();

  static int ng_dens_hse;
  static int ng_pres_hse;

  //! Initialize HSE
  void initHSE();
  void erf_enforce_hse(amrex::Vector<amrex::Real>& dens,
                       amrex::Vector<amrex::Real>& pres);

  //! Initialize grid data at problem start-up.
  virtual void initData() override;
  void initDataProb(amrex::MultiFab& S_new,
                    amrex::MultiFab& U_new,
                    amrex::MultiFab& V_new,
                    amrex::MultiFab& W_new);

  //! Get problem metrics.
  amrex::MultiFab* Area();
  amrex::MultiFab& Area(int dir);

  //! The volume of the multifab.
  amrex::MultiFab& Volume();

  //! Set time levels of state data.
  virtual void setTimeLevel(
    amrex::Real time, amrex::Real dt_old, amrex::Real dt_new) override;

  //! Initialize data on this level from another ERF (during regrid).
  virtual void init(amrex::AmrLevel& old) override;

  /** Initialize data on this level after regridding if old level did not
   *  previously exist
   */
  virtual void init() override;

  //! for the Amr class to do timed load balances.
  virtual int WorkEstType() override { return Work_Estimate_Type; }

  //! Proceed with next timestep?
  virtual int okToContinue() override;

  //! Advance grids at this level in time.
  virtual amrex::Real
  advance(amrex::Real time, amrex::Real dt, int iteration, int ncycle) override;

  //! Estimate time step.
  amrex::Real estTimeStep(amrex::Real dt_old);

  //! Compute initial time step.
  amrex::Real initialTimeStep();

  //! Compute initial `dt'.
  virtual void computeInitialDt(
    int finest_level,
    int sub_cycle,
    amrex::Vector<int>& n_cycle,
    const amrex::Vector<amrex::IntVect>& ref_ratio,
    amrex::Vector<amrex::Real>& dt_level,
    amrex::Real stop_time) override;

  //! Compute new `dt'.
  virtual void computeNewDt(
    int finest_level,
    int sub_cycle,
    amrex::Vector<int>& n_cycle,
    const amrex::Vector<amrex::IntVect>& ref_ratio,
    amrex::Vector<amrex::Real>& dt_min,
    amrex::Vector<amrex::Real>& dt_level,
    amrex::Real stop_time,
    int post_regrid_flag) override;

  //! Allocate data at old time.
  virtual void allocOldData() override;

  //! Remove data at old time.
  virtual void removeOldData() override;

  //! Passes some data about the grid to a Fortran module.
  void setGridInfo();

  // Print information about energy budget.
  //void do_energy_diagnostics();

  //! Do work after timestep().
  virtual void post_timestep(int iteration) override;

  //! Contains operations to be done only after a full coarse timestep.
  virtual void postCoarseTimeStep(amrex::Real cumtime) override;

  //! Do work after regrid().
  virtual void post_regrid(int lbase, int new_finest) override;

  //! Do work after a restart().
  virtual void post_restart() override;

  //! Do work after init().
  virtual void post_init(amrex::Real stop_time) override;

  //! Error estimation for regridding.
  virtual void errorEst(
    amrex::TagBoxArray& tb,
    int clearval,
    int tagval,
    amrex::Real time,
    int n_error_buf = 0,
    int ngrow = 0) override;

  /** Returns a MultiFab containing the derived data for this level.
   * The user is responsible for deleting this pointer when done
   * with it.  If ngrow>0 the MultiFab is built on the appropriately
   * grown BoxArray.
   */
  virtual std::unique_ptr<amrex::MultiFab>
  derive(const std::string& name, amrex::Real time, int ngrow) override;

  /** This version of derive() fills the dcomp'th component of mf with the
   * derived quantity.
   */
  virtual void derive(
    const std::string& name,
    amrex::Real time,
    amrex::MultiFab& mf,
    int dcomp) override;

  static int numGrow();

  static void applyBCs( const amrex::Geometry geom, amrex::Vector<MultiFab*>& vars) {
    for (int dir = 0; dir < 2*AMREX_SPACEDIM; dir++) {
      bc_recs[dir]->applyBC(geom, vars);
    }
  }

  amrex::Real volWgtSum( const std::string& name,
                         amrex::Real time,
                         bool local = false,
                         bool finemask = true);

  amrex::Real volWgtSquaredSum( const std::string& name, amrex::Real time, bool local = false);

  amrex::Real volWgtSumMF( const amrex::MultiFab& mf,
                           int comp,
                           bool local = false,
                           bool finemask = false);

  amrex::Real volWgtSquaredSumDiff( int comp, amrex::Real time, bool local = false);
  amrex::Real maxDerive( const std::string& name, amrex::Real time, bool local = false);

  static int FirstAdv, NumAdv;

  // This MultiFab is on the coarser level.  This is useful for the coarser
  // level to mask out the finer level.  We only build this when it is needed.
  // This coarse MultiFab has to live on the fine level because it must be
  // updated even when only the fine level changes.
  amrex::MultiFab fine_mask;
  amrex::MultiFab& build_fine_mask();

  amrex::YAFluxRegister flux_reg;

  amrex::YAFluxRegister& getFluxReg();
  amrex::YAFluxRegister& getFluxReg(int lev);

  void avgDown();
  void avgDown(int state_indx);

protected:
  friend class IOManager;

  // io manager
  std::unique_ptr<IOManager> io_mgr;

  // Build a mask that ghost cells overlapping with interior cells in the same
  // multifab are set to 0, whereas others are set to 1.
  amrex::Vector<std::unique_ptr<amrex::iMultiFab>> ib_mask;
  const amrex::iMultiFab* build_interior_boundary_mask(int ng);

  /// A state array with ghost zones.
  amrex::MultiFab Sborder;

  static void clear_prob();

  static void read_params();

  static void read_tagging_params();

  template<int DIM, math_bcs::BCBound Bound>
  static std::unique_ptr<phys_bcs::BCBase> initialize_bcs(const std::string& bc_char);

  ERF& getLevel(int lev);

  void reflux();

  void buildMetrics();

  /// integrate derived quantities over domain
  amrex::Real
  sumDerive(const std::string& name, amrex::Real time, bool local = false);

  void sum_integrated_quantities();

  //void write_info();

  void stopJob();

  // ***************************************************************************
  // The data or the member variables
  amrex::MultiFab volume;
  amrex::MultiFab area[AMREX_SPACEDIM];
  int m_ngrow;

  // ***************************************************************************
  // Static data members.
  static bool signalStopJob;
  static bool dump_old;
  static int verbose;

  static std::string coupling_type;
  static int do_reflux;
  static int do_avg_down;

  static int sum_interval;
  static amrex::Real sum_per;

  static amrex::Real initial_dt;
  static amrex::Real fixed_dt;
  static amrex::Real max_dt;
  static amrex::Real dt_cutoff;
  static amrex::Real init_shrink;
  static amrex::Real change_max;
  static amrex::Real cfl;
  static amrex::Vector<std::unique_ptr<phys_bcs::BCBase> > bc_recs;

  /** for keeping track of the amount of CPU time used -- this will persist
   * after restarts
   */
  static amrex::Real previousCPUTimeUsed;
  static amrex::Real startCPUTime;

  static amrex::Real getCPUTime();

  bool FillPatchedOldState_ok;

  /// There can be only one Diffusion object, it covers all levels.
  static class Diffusion* diffusion;

  static int num_state_type;

  static amrex::Vector<int> src_list;

  // Base state quantities: density and pressure, 1 array per AMR level
  static amrex::Vector<amrex::Vector<amrex::Real> > h_dens_hse;
  static amrex::Vector<amrex::Vector<amrex::Real> > h_pres_hse;
  static amrex::Vector<amrex::Gpu::DeviceVector<amrex::Real> > d_dens_hse;
  static amrex::Vector<amrex::Gpu::DeviceVector<amrex::Real> > d_pres_hse;

/* problem-specific includes */
#include <Problem.H>

};

void erf_bcfill_hyp(
  amrex::Box const& bx,
  amrex::FArrayBox& data,
  const int dcomp,
  const int numcomp,
  amrex::Geometry const& geom,
  const amrex::Real time,
  const amrex::Vector<amrex::BCRec>& bcr,
  const int bcomp,
  const int scomp);

void erf_nullfill(
  amrex::Box const& bx,
  amrex::FArrayBox& data,
  const int dcomp,
  const int numcomp,
  amrex::Geometry const& geom,
  const amrex::Real time,
  const amrex::Vector<amrex::BCRec>& bcr,
  const int bcomp,
  const int scomp);

// *****************************************************************************
// Inlines.
//

AMREX_FORCE_INLINE
int
ERF::numGrow()
{
  return NUM_GROW;
}

AMREX_FORCE_INLINE
amrex::MultiFab*
ERF::Area()
{
  return area;
}

AMREX_FORCE_INLINE
amrex::MultiFab&
ERF::Area(int dir)
{
  return area[dir];
}

AMREX_FORCE_INLINE
amrex::MultiFab&
ERF::Volume()
{
  return volume;
}

AMREX_FORCE_INLINE
ERF&
ERF::getLevel(int lev)
{
  return *(ERF*)&parent->getLevel(lev);
}

AMREX_FORCE_INLINE
void
ERF::stopJob()
{
  signalStopJob = true;
}

AMREX_FORCE_INLINE
amrex::YAFluxRegister&
ERF::getFluxReg()
{
  return flux_reg;
}

AMREX_FORCE_INLINE
amrex::YAFluxRegister&
ERF::getFluxReg(int lev)
{
  return getLevel(lev).getFluxReg();
}
#endif
