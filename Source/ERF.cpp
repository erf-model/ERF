#ifdef _OPENMP
#include <omp.h>
#endif

#include <AMReX_Vector.H>
#include <AMReX_Utility.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_VisMF.H>
#include <AMReX_TagBox.H>
#include <AMReX_ParmParse.H>

#ifdef ERF_USE_MASA
#include <masa.h>
using namespace MASA;
#endif

#include "ERF.H"
#include "Derive.H"
#include "prob.H"
#include "Transport.H"
#include "EOS.H"
#include "Utilities.H"
#include "Tagging.H"
#include "IndexDefines.H"

using namespace amrex;

bool ERF::signalStopJob = false;
bool ERF::dump_old = false;
int ERF::verbose = 0;
amrex::BCRec ERF::phys_bc;
amrex::Real ERF::frac_change = 1.e200;
int ERF::NumAdv = 0;
int ERF::FirstAdv = -1;

#include "erf_defaults.H"

int ERF::nGrowTr = 4;

bool ERF::do_mol_load_balance = false;

amrex::Vector<int> ERF::src_list;

// this will be reset upon restart
amrex::Real ERF::previousCPUTimeUsed = 0.0;
amrex::Real ERF::startCPUTime = 0.0;
int ERF::num_state_type = 0;

void
ERF::variableCleanUp()
{
  desc_lst.clear();

  transport_close();

  clear_prob();
}

void
ERF::read_params()
{
  static bool read_params_done = false;

  if (read_params_done)
    return;

  read_params_done = true;

  amrex::ParmParse pp("erf");

#include <erf_queries.H>

  pp.query("v", verbose);
  pp.query("sum_interval", sum_interval);
  pp.query("dump_old", dump_old);

  // Get boundary conditions
  amrex::Vector<std::string> lo_bc_char(AMREX_SPACEDIM);
  amrex::Vector<std::string> hi_bc_char(AMREX_SPACEDIM);
  pp.getarr("lo_bc", lo_bc_char, 0, AMREX_SPACEDIM);
  pp.getarr("hi_bc", hi_bc_char, 0, AMREX_SPACEDIM);

  amrex::Vector<int> lo_bc(AMREX_SPACEDIM), hi_bc(AMREX_SPACEDIM);
  for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
    if (!lo_bc_char[dir].compare("Interior")) {
      lo_bc[dir] = 0;
      std::cout << "DIR IS INTERIOR " << dir << std::endl;
    } else if (!lo_bc_char[dir].compare("Hard")) {
      lo_bc[dir] = 1;
    } else if (!lo_bc_char[dir].compare("FOExtrap")) {
      lo_bc[dir] = 2;
    } else if (!lo_bc_char[dir].compare("Symmetry")) {
      lo_bc[dir] = 3;
    } else if (!lo_bc_char[dir].compare("SlipWall")) {
      lo_bc[dir] = 4;
    } else if (!lo_bc_char[dir].compare("NoSlipWall")) {
      lo_bc[dir] = 5;
    } else if (!lo_bc_char[dir].compare("UserBC")) {
      lo_bc[dir] = 6;
    } else {
      amrex::Abort("Wrong boundary condition word in lo_bc, please use: "
                   "Interior, UserBC, Symmetry, SlipWall, NoSlipWall");
    }

    if (!hi_bc_char[dir].compare("Interior")) {
      hi_bc[dir] = 0;
    } else if (!hi_bc_char[dir].compare("Hard")) {
      hi_bc[dir] = 1;
    } else if (!hi_bc_char[dir].compare("FOExtrap")) {
      hi_bc[dir] = 2;
    } else if (!hi_bc_char[dir].compare("Symmetry")) {
      hi_bc[dir] = 3;
    } else if (!hi_bc_char[dir].compare("SlipWall")) {
      hi_bc[dir] = 4;
    } else if (!hi_bc_char[dir].compare("NoSlipWall")) {
      hi_bc[dir] = 5;
    } else if (!hi_bc_char[dir].compare("UserBC")) {
      hi_bc[dir] = 6;
    } else {
      amrex::Abort("Wrong boundary condition word in hi_bc, please use: "
                   "Interior, UserBC, Symmetry, SlipWall, NoSlipWall");
    }
  }

  for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
    phys_bc.setLo(dir, lo_bc[dir]);
    phys_bc.setHi(dir, hi_bc[dir]);
  }

  //
  // Check phys_bc against possible periodic geometry
  // if periodic, must have internal BC marked.
  //
  //
  // Do idiot check.  Periodic means interior in those directions.
  //
  for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
    if (amrex::DefaultGeometry().isPeriodic(dir)) {
      if (
        lo_bc[dir] != Interior && amrex::ParallelDescriptor::IOProcessor()) {
        std::cerr << "ERF::read_params:periodic in direction " << dir
                  << " but low BC is not Interior\n";
        amrex::Error();
      }
      if (
        hi_bc[dir] != Interior && amrex::ParallelDescriptor::IOProcessor()) {
        std::cerr << "ERF::read_params:periodic in direction " << dir
                  << " but high BC is not Interior\n";
        amrex::Error();
      }
    } else {
      //
      // Do idiot check. If not periodic, should not be interior.
      //
      if (lo_bc[dir] == Interior && amrex::ParallelDescriptor::IOProcessor()) {
        std::cerr << "ERF::read_params:interior bc in direction " << dir
                  << " but not periodic\n";
        amrex::Error();
      }
      if (hi_bc[dir] == Interior && amrex::ParallelDescriptor::IOProcessor()) {
        std::cerr << "ERF::read_params:interior bc in direction " << dir
                  << " but not periodic\n";
        amrex::Error();
      }
    }
  }

  if (amrex::DefaultGeometry().IsRZ() && (lo_bc[0] != Symmetry)) {
    amrex::Error("ERF::read_params: must set r=0 boundary condition to "
                 "Symmetry for r-z");
  }

  // TODO: Any reason to support spherical in ERF?
  if (amrex::DefaultGeometry().IsRZ()) {
    amrex::Abort("We don't support cylindrical coordinate systems in 3D");
  } else if (amrex::DefaultGeometry().IsSPHERICAL()) {
    amrex::Abort("We don't support spherical coordinate systems in 3D");
  }

  // sanity checks
  if (cfl <= 0.0 || cfl > 1.0) {
    amrex::Error("Invalid CFL factor; must be between zero and one.");
  }

  if (max_dt < fixed_dt) {
    amrex::Error("Cannot have max_dt < fixed_dt");
  }

  // Read tagging parameters
  read_tagging_params();

  // TODO: What is this?
  amrex::StateDescriptor::setBndryFuncThreadSafety(bndry_func_thread_safe);

  // Get some useful amr inputs
  amrex::ParmParse ppa("amr");

  // This turns on the lb stuff inside Amr, but we use our own flag to signal
  // whether to gather data
  ppa.query("loadbalance_with_workestimates", do_mol_load_balance);
}

ERF::ERF()
  : old_sources(num_src),
    new_sources(num_src)
#ifdef ERF_USE_MASA
    ,
    mms_src_evaluated(false)
#endif
{
}

ERF::ERF(
  amrex::Amr& papa,
  int lev,
  const amrex::Geometry& level_geom,
  const amrex::BoxArray& bl,
  const amrex::DistributionMapping& dm,
  amrex::Real time)
  : AmrLevel(papa, lev, level_geom, bl, dm, time),
    old_sources(num_src),
    new_sources(num_src)
#ifdef ERF_USE_MASA
    ,
    mms_src_evaluated(false)
#endif
{
  buildMetrics();

  amrex::MultiFab& S_new = get_new_data(State_Type);

  for (int n = 0; n < src_list.size(); ++n) {
    int oldGrow = NUM_GROW;
    int newGrow = S_new.nGrow();
    old_sources[src_list[n]] =
      std::unique_ptr<amrex::MultiFab>(new amrex::MultiFab(
        grids, dmap, NVAR, oldGrow, amrex::MFInfo(), Factory()));
    new_sources[src_list[n]] =
      std::unique_ptr<amrex::MultiFab>(new amrex::MultiFab(
        grids, dmap, NVAR, newGrow, amrex::MFInfo(), Factory()));
  }

  Sborder.define(grids, dmap, NVAR, NUM_GROW, amrex::MFInfo(), Factory());

  if (do_reflux && level > 0) {
    flux_reg.define(
      bl, papa.boxArray(level - 1), dm, papa.DistributionMap(level - 1),
      level_geom, papa.Geom(level - 1), papa.refRatio(level - 1), level, NVAR);
  }
}

ERF::~ERF() {}

void
ERF::buildMetrics()
{
  const int ngrd = grids.size();

  const amrex::Real* dx = geom.CellSize();

  volume.clear();
  volume.define(
    grids, dmap, 1, NUM_GROW, amrex::MFInfo(), amrex::FArrayBoxFactory());
  geom.GetVolume(volume);

  for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
    area[dir].clear();
    area[dir].define(
      getEdgeBoxArray(dir), dmap, 1, NUM_GROW, amrex::MFInfo(),
      amrex::FArrayBoxFactory());
    geom.GetFaceArea(area[dir], dir);
  }

  level_mask.clear();
  level_mask.define(grids, dmap, 1, 1);
  level_mask.BuildMask(
    geom.Domain(), geom.periodicity(), levmsk_covered, levmsk_notcovered,
    levmsk_physbnd, levmsk_interior);

  if (level == 0)
    setGridInfo();
}

void
ERF::setTimeLevel(amrex::Real time, amrex::Real dt_old, amrex::Real dt_new)
{
  AmrLevel::setTimeLevel(time, dt_old, dt_new);
}

void
ERF::setGridInfo()
{
  /** Send refinement data to Fortran. We do it here
    because now the grids have been initialized and
    we need this data for setting up the problem.
    Note that this routine will always get called
    on level 0, even if we are doing a restart,
    so it is safe to put this here.
    */

  if (level == 0) {
    int max_level = parent->maxLevel();
    int nlevs = max_level + 1;

    amrex::Real dx_level[3 * nlevs];
    int domlo_level[3 * nlevs];
    int domhi_level[3 * nlevs];

    const amrex::Real* dx_coarse = geom.CellSize();

    const int* domlo_coarse = geom.Domain().loVect();
    const int* domhi_coarse = geom.Domain().hiVect();

    for (int dir = 0; dir < 3; dir++) {
      dx_level[dir] = (ZFILL(dx_coarse))[dir];

      domlo_level[dir] = (ARLIM_3D(domlo_coarse))[dir];
      domhi_level[dir] = (ARLIM_3D(domhi_coarse))[dir];
    }

    for (int lev = 1; lev <= max_level; lev++) {
      amrex::IntVect ref_ratio = parent->refRatio(lev - 1);

      // Note that we are explicitly calculating here what the
      // data would be on refined levels rather than getting the
      // data directly from those levels, because some potential
      // refined levels may not exist at the beginning of the simulation.

      for (int dir = 0; dir < 3; dir++) {
        if (dir < AMREX_SPACEDIM) {
          dx_level[3 * lev + dir] =
            dx_level[3 * (lev - 1) + dir] / ref_ratio[dir];
          int ncell = (domhi_level[3 * (lev - 1) + dir] -
                       domlo_level[3 * (lev - 1) + dir] + 1) *
                      ref_ratio[dir];
          domlo_level[3 * lev + dir] = domlo_level[dir];
          domhi_level[3 * lev + dir] = domlo_level[3 * lev + dir] + ncell - 1;
        } else {
          dx_level[3 * lev + dir] = 0.0;
          domlo_level[3 * lev + dir] = 0;
          domhi_level[3 * lev + dir] = 0;
        }
      }
    }

    // Don't need this in pure C++?
    // set_grid_info(max_level, dx_level, domlo_level, domhi_level);
  }
}

/*
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
erf_prob_close()
{
}
*/

void
ERF::initData()
{
  BL_PROFILE("ERF::initData()");

  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();
  amrex::Real cur_time = state[State_Type].curTime();

  amrex::MultiFab& S_new = get_new_data(State_Type);
  amrex::MultiFab& U_new = get_new_data(X_Vel_Type);
  amrex::MultiFab& V_new = get_new_data(Y_Vel_Type);
  amrex::MultiFab& W_new = get_new_data(Z_Vel_Type);

  // Initialize to zero (though we sholdn't actually need to do this)
  S_new.setVal(0.0);
  U_new.setVal(0.0);
  V_new.setVal(0.0);
  W_new.setVal(0.0);

  if (verbose) {
    amrex::Print() << "Initializing the data at level " << level << std::endl;
  }

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
  for (amrex::MFIter mfi(S_new, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) 
  {
    const amrex::Box& bx = mfi.tilebox();
    auto sfab  = S_new.array(mfi);
    auto ufab  = U_new.array(mfi);
    auto vfab  = V_new.array(mfi);
    auto wfab  = W_new.array(mfi);
    const auto geomdata = geom.data();

    // Construct a box that is on x-faces
    const amrex::Box& xbx = surroundingNodes(bx,0);

    // Call for all (i,j,k) in the x-face-centered box
    amrex::ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      erf_init_xvel(i, j, k, ufab, geomdata);
    });

    // Construct a box that is on y-faces
    const amrex::Box& ybx = surroundingNodes(bx,1);

    // Call for all (i,j,k) in the y-face-centered box
    amrex::ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      erf_init_yvel(i, j, k, vfab, geomdata);
    });

    // Construct a box that is on z-faces
    const amrex::Box& zbx = surroundingNodes(bx,2);

    // Call for all (i,j,k) in the z-face-centered box
    amrex::ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      erf_init_zvel(i, j, k, wfab, geomdata);
    });
  }

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
  for (amrex::MFIter mfi(S_new, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) 
  {
    const amrex::Box& bx = mfi.tilebox();
    auto sfab  = S_new.array(mfi);
    auto ufab  = U_new.array(mfi);
    auto vfab  = V_new.array(mfi);
    auto wfab  = W_new.array(mfi);
    const auto geomdata = geom.data();

    // Call for all (i,j,k) in the cell-centered box
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      erf_initdata(i, j, k, sfab, ufab, vfab, wfab, geomdata);
    });
  }

  if (verbose) {
    amrex::Print() << "Done initializing level " << level << " data "
                   << std::endl;
  }
}

void
ERF::init(AmrLevel& old)
{
  BL_PROFILE("ERF::init(old)");

  ERF* oldlev = (ERF*)&old;

  //
  // Create new grid data by fillpatching from old.
  //
  amrex::Real dt_new = parent->dtLevel(level);
  amrex::Real cur_time = oldlev->state[State_Type].curTime();
  amrex::Real prev_time = oldlev->state[State_Type].prevTime();
  amrex::Real dt_old = cur_time - prev_time;
  setTimeLevel(cur_time, dt_old, dt_new);

  amrex::MultiFab& S_new = get_new_data(State_Type);
  FillPatch(old, S_new, 0, cur_time, State_Type, 0, NVAR);

  if (do_mol_load_balance) {
    amrex::MultiFab& work_estimate_new = get_new_data(Work_Estimate_Type);
    FillPatch(
      old, work_estimate_new, 0, cur_time, Work_Estimate_Type, 0,
      work_estimate_new.nComp());
  }
}

void
ERF::init()
{
  /**
     This version inits the data on a new level that did not
     exist before regridding.
  */
  BL_PROFILE("ERF::init()");

  amrex::Real dt = parent->dtLevel(level);
  amrex::Real cur_time = getLevel(level - 1).state[State_Type].curTime();
  amrex::Real prev_time = getLevel(level - 1).state[State_Type].prevTime();

  amrex::Real dt_old =
    (cur_time - prev_time) / (amrex::Real)parent->MaxRefRatio(level - 1);

  setTimeLevel(cur_time, dt_old, dt);
  amrex::MultiFab& S_new = get_new_data(State_Type);
  FillCoarsePatch(S_new, 0, cur_time, State_Type, 0, NVAR);

  if (do_mol_load_balance) {
    amrex::MultiFab& work_estimate_new = get_new_data(Work_Estimate_Type);
    int ncomp = work_estimate_new.nComp();
    FillCoarsePatch(
      work_estimate_new, 0, cur_time, Work_Estimate_Type, 0, ncomp);
  }
}

amrex::Real
ERF::initialTimeStep()
{
  BL_PROFILE("ERF::initialTimeStep()");

  amrex::Real dummy_dt = 0.0;
  amrex::Real init_dt = 0.0;

  if (initial_dt > 0.0) {
    init_dt = initial_dt;
  } else {
    init_dt = init_shrink * estTimeStep(dummy_dt);
  }

  return init_dt;
}

amrex::Real
ERF::estTimeStep(amrex::Real dt_old)
{
  BL_PROFILE("ERF::estTimeStep()");

  if (fixed_dt > 0.0)
    return fixed_dt;

  amrex::Real estdt = max_dt;

  const amrex::Real* dx = geom.CellSize();

  std::string limiter = "erf.max_dt";

  // Start the hydro with the max_dt value, but divide by CFL
  // to account for the fact that we multiply by it at the end.
  // This ensures that if max_dt is more restrictive than the hydro
  // criterion, we will get exactly max_dt for a timestep.

  amrex::Real estdt_hydro = max_dt / cfl;

  auto const dxinv = geom.InvCellSizeArray();

  MultiFab& state = get_new_data(State_Type);

  MultiFab const& x_vel_on_face = get_new_data(X_Vel_Type);
  MultiFab const& y_vel_on_face = get_new_data(Y_Vel_Type);
  MultiFab const& z_vel_on_face = get_new_data(Z_Vel_Type);

  MultiFab ccvel(grids,dmap,3,0);

  average_face_to_cellcenter(ccvel,0,
      Array<const MultiFab*,3>{&x_vel_on_face, &y_vel_on_face, &z_vel_on_face});

  // Initialize to large value since we are taking min below
  // Real estdt_hydro_inv = 1.e100;

  Real estdt_hydro_inv = amrex::ReduceMax(state, ccvel, 0,
       [=] AMREX_GPU_HOST_DEVICE (Box const& b,
                                  Array4<Real const> const& s,
                                  Array4<Real const> const& u) -> Real
       {
           Real dt = -1.e100;
           amrex::Loop(b, [=,&dt] (int i, int j, int k) noexcept
           {
               const amrex::Real rho   = s(i, j, k, Density_comp);
               const amrex::Real theta = s(i, j, k,   Theta_comp);

               amrex::Real pressure = getPgivenRTh(rho,theta);
               amrex::Real c = std::sqrt(Gamma * pressure / rho);

               dt = amrex::max(((amrex::Math::abs(u(i,j,k,0))+c)*dxinv[0]),
                               ((amrex::Math::abs(u(i,j,k,1))+c)*dxinv[1]),
                               ((amrex::Math::abs(u(i,j,k,2))+c)*dxinv[2]), dt);
           });
           return dt;
       });

   amrex::ParallelDescriptor::ReduceRealMax(estdt_hydro_inv);
   estdt_hydro = cfl / estdt_hydro_inv;;

   if (verbose) {
     amrex::Print() << "...estimated hydro-limited timestep at level " << level
                    << ": " << estdt_hydro << std::endl;
   }

    // Determine if this is more restrictive than the maximum timestep limiting
    if (estdt_hydro < estdt) {
      limiter = "hydro";
      estdt = estdt_hydro;
    }

  if (verbose) {
    amrex::Print() << "ERF::estTimeStep (" << limiter << "-limited) at level "
                   << level << ":  estdt = " << estdt << '\n';
  }

  return estdt;
}

void
ERF::computeNewDt(
  int finest_level,
  int sub_cycle,
  amrex::Vector<int>& n_cycle,
  const amrex::Vector<amrex::IntVect>& ref_ratio,
  amrex::Vector<amrex::Real>& dt_min,
  amrex::Vector<amrex::Real>& dt_level,
  amrex::Real stop_time,
  int post_regrid_flag)
{
  BL_PROFILE("ERF::computeNewDt()");

  //
  // We are at the start of a coarse grid timecycle.
  // Compute the timesteps for the next iteration.
  //
  if (level > 0)
    return;

  amrex::Real dt_0 = 1.0e+100;
  int n_factor = 1;
  for (int i = 0; i <= finest_level; i++) {
    ERF& adv_level = getLevel(i);
    dt_min[i] = adv_level.estTimeStep(dt_level[i]);
  }

  if (fixed_dt <= 0.0) {
    if (post_regrid_flag == 1) {
      // Limit dt's by pre-regrid dt
      for (int i = 0; i <= finest_level; i++) {
        dt_min[i] = amrex::min(dt_min[i], dt_level[i]);
      }
    } else {
      // Limit dt's by change_max * old dt
      for (int i = 0; i <= finest_level; i++) {
        if (verbose && amrex::ParallelDescriptor::IOProcessor()) {
          if (dt_min[i] > change_max * dt_level[i]) {
            amrex::Print() << "ERF::compute_new_dt : limiting dt at level "
                           << i << '\n';
            amrex::Print() << " ... new dt computed: " << dt_min[i] << '\n';
            amrex::Print() << " ... but limiting to: "
                           << change_max * dt_level[i] << " = " << change_max
                           << " * " << dt_level[i] << '\n';
          }
        }
        dt_min[i] = amrex::min(dt_min[i], change_max * dt_level[i]);
      }
    }
  }

  // Find the minimum over all levels
  for (int i = 0; i <= finest_level; i++) {
    n_factor *= n_cycle[i];
    dt_0 = amrex::min(dt_0, n_factor * dt_min[i]);
  }

  // Limit dt's by the value of stop_time.
  const amrex::Real dt_eps = 0.001 * dt_0;
  amrex::Real cur_time = state[State_Type].curTime();
  if (stop_time >= 0.0) {
    if ((cur_time + dt_0) > (stop_time - dt_eps)) {
      dt_0 = stop_time - cur_time;
    }
  }

  n_factor = 1;
  for (int i = 0; i <= finest_level; i++) {
    n_factor *= n_cycle[i];
    dt_level[i] = dt_0 / n_factor;
  }
}

void
ERF::computeInitialDt(
  int finest_level,
  int sub_cycle,
  amrex::Vector<int>& n_cycle,
  const amrex::Vector<amrex::IntVect>& ref_ratio,
  amrex::Vector<amrex::Real>& dt_level,
  amrex::Real stop_time)
{
  BL_PROFILE("ERF::computeInitialDt()");

  // Grids have been constructed, compute dt for all levels.
  if (level > 0)
    return;

  amrex::Real dt_0 = 1.0e+100;
  int n_factor = 1;
  /// TODO/DEBUG: This will need to change for optimal subcycling.
  for (int i = 0; i <= finest_level; i++) {
    dt_level[i] = getLevel(i).initialTimeStep();
    n_factor *= n_cycle[i];
    dt_0 = amrex::min(dt_0, n_factor * dt_level[i]);
  }

  // Limit dt's by the value of stop_time.
  const amrex::Real dt_eps = 0.001 * dt_0;
  amrex::Real cur_time = state[State_Type].curTime();
  if (stop_time >= 0.0) {
    if ((cur_time + dt_0) > (stop_time - dt_eps)) {
      dt_0 = stop_time - cur_time;
    }
  }

  n_factor = 1;
  for (int i = 0; i <= finest_level; i++) {
    n_factor *= n_cycle[i];
    dt_level[i] = dt_0 / n_factor;
  }
}

void
ERF::post_timestep(int iteration)
{
  BL_PROFILE("ERF::post_timestep()");

  const int finest_level = parent->finestLevel();
  const int ncycle = parent->nCycle(level);

  if (do_reflux && level < finest_level) {
    reflux();

    // We need to do this before anything else because refluxing changes the
    // values of coarse cells
    //    underneath fine grids with the assumption they'll be over-written by
    //    averaging down
    if (level < finest_level) {
      avgDown();
    }

    // Clean up any aberrant state data generated by the reflux.
    amrex::MultiFab& S_new_crse = get_new_data(State_Type);
    // clean_state(S_new_crse);
  }

  // Re-compute temperature after all the other updates.
  amrex::MultiFab& S_new = get_new_data(State_Type);
  int ng_pts = 0;

#ifdef DO_PROBLEM_POST_TIMESTEP

  problem_post_timestep();

#endif

  if (level == 0) {
    int nstep = parent->levelSteps(0);
    amrex::Real dtlev = parent->dtLevel(0);
    amrex::Real cumtime = parent->cumTime() + dtlev;

    bool sum_int_test = (sum_interval > 0 && nstep % sum_interval == 0);

    bool sum_per_test = false;

    if (sum_per > 0.0) {
      const int num_per_old = amrex::Math::floor((cumtime - dtlev) / sum_per);
      const int num_per_new = amrex::Math::floor((cumtime) / sum_per);

      if (num_per_old != num_per_new) {
        sum_per_test = true;
      }
    }

    if (sum_int_test || sum_per_test) {
      sum_integrated_quantities();
    }
  }
}

void
ERF::post_restart()
{
  BL_PROFILE("ERF::post_restart()");

  amrex::Real cur_time = state[State_Type].curTime();

#ifdef DO_PROBLEM_POST_RESTART
  problem_post_restart();
#endif
}

void
ERF::postCoarseTimeStep(amrex::Real cumtime)
{
  BL_PROFILE("ERF::postCoarseTimeStep()");
  AmrLevel::postCoarseTimeStep(cumtime);
}

void
ERF::post_regrid(int lbase, int new_finest)
{
  BL_PROFILE("ERF::post_regrid()");
  fine_mask.clear();
}

void
ERF::post_init(amrex::Real stop_time)
{
  BL_PROFILE("ERF::post_init()");

  amrex::Real dtlev = parent->dtLevel(level);
  amrex::Real cumtime = parent->cumTime();

  if (level > 0)
    return;

  //
  // Average data down from finer levels
  // so that conserved data is consistent between levels.
  //
  if (do_avg_down != 0) {
    int finest_level = parent->finestLevel();
    for (int k = finest_level - 1; k >= 0; k--) {
      getLevel(k).avgDown();
    }
  }

#ifdef DO_PROBLEM_POST_INIT
  //
  // Allow the user to define their own post_init functions.
  //
  problem_post_init();
#endif

  int nstep = parent->levelSteps(0);
  if (cumtime != 0.0)
    cumtime += dtlev;

  bool sum_int_test = false;

  if (sum_interval > 0) {
    if (nstep % sum_interval == 0) {
      sum_int_test = true;
    }
  }

  bool sum_per_test = false;

  if (sum_per > 0.0) {
    const int num_per_old = amrex::Math::floor((cumtime - dtlev) / sum_per);
    const int num_per_new = amrex::Math::floor((cumtime) / sum_per);

    if (num_per_old != num_per_new) {
      sum_per_test = true;
    }
  }

  if (sum_int_test || sum_per_test) {
    sum_integrated_quantities();
  }
}

int
ERF::okToContinue()
{
  if (level > 0) {
    return 1;
  }

  int test = 1;

  if (signalStopJob) {
    test = 0;

    amrex::Print()
      << " Signalling a stop of the run due to signalStopJob = true."
      << std::endl;
  } else if (parent->dtLevel(0) < dt_cutoff) {
    test = 0;

    amrex::Print() << " Signalling a stop of the run because dt < dt_cutoff."
                   << std::endl;
  }

  return test;
}

void
ERF::reflux()
{
  BL_PROFILE("ERF::reflux()");

  AMREX_ASSERT(level < parent->finestLevel());

  const amrex::Real strt = amrex::ParallelDescriptor::second();

  ERF& fine_level = getLevel(level + 1);
  amrex::MultiFab& S_crse = get_new_data(State_Type);

  fine_level.flux_reg.Reflux(S_crse);

  if (verbose) {
    const int IOProc = amrex::ParallelDescriptor::IOProcessorNumber();
    amrex::Real end = amrex::ParallelDescriptor::second() - strt;

#ifdef AMREX_LAZY
    Lazy::QueueReduction([=]() mutable {
#endif
      amrex::ParallelDescriptor::ReduceRealMax(end, IOProc);

      amrex::Print() << "ERF::reflux() at level " << level
                     << " : time = " << end << std::endl;
#ifdef AMREX_LAZY
    });
#endif
  }
}

void
ERF::avgDown()
{
  BL_PROFILE("ERF::avgDown()");

  if (level == parent->finestLevel())
    return;

  avgDown(State_Type);
}

void
ERF::avgDown(int state_indx)
{
  BL_PROFILE("ERF::avgDown(state_indx)");

  if (level == parent->finestLevel())
    return;

  amrex::MultiFab& S_crse = get_new_data(state_indx);
  amrex::MultiFab& S_fine = getLevel(level + 1).get_new_data(state_indx);

  const amrex::Geometry& fgeom = getLevel(level + 1).geom;
  const amrex::Geometry& cgeom = geom;

  amrex::average_down(
    S_fine, S_crse, fgeom, cgeom, 0, S_fine.nComp(), fine_ratio);
}

void
ERF::allocOldData()
{
  for (int k = 0; k < num_state_type; k++) {
    state[k].allocOldData();
  }
}

void
ERF::removeOldData()
{
  AmrLevel::removeOldData();
}

void
ERF::errorEst(
  amrex::TagBoxArray& tags,
  int /*clearval*/,
  int /*tagval*/,
  amrex::Real time,
  int n_error_buf,
  int ngrow)
{
  BL_PROFILE("ERF::errorEst()");

  amrex::MultiFab S_data(
    get_new_data(State_Type).boxArray(),
    get_new_data(State_Type).DistributionMap(), NVAR, 1);
  const amrex::Real cur_time = state[State_Type].curTime();
  FillPatch(
    *this, S_data, S_data.nGrow(), cur_time, State_Type, Density_comp, NVAR, 0);

  amrex::Vector<amrex::BCRec> bcs(NVAR);
  const char tagval = amrex::TagBox::SET;
  const char clearval = amrex::TagBox::CLEAR;

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
  {
    for (amrex::MFIter mfi(S_data, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {
      const amrex::Box& tilebox = mfi.tilebox();
      const auto Sfab = S_data.array(mfi);
      auto tag_arr = tags.array(mfi);
      const auto datbox = S_data[mfi].box();
      amrex::Elixir S_data_mfi_eli = S_data[mfi].elixir();

      amrex::FArrayBox S_derData(datbox, 1);
      amrex::Elixir S_derData_eli = S_derData.elixir();
      auto S_derarr = S_derData.array();
      const int ncp = S_derData.nComp();
      const int* bc = bcs[0].data();

      // Tagging density
      if (level < TaggingParm::max_denerr_lev) {
        amrex::ParallelFor(
          tilebox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            tag_error(i, j, k, tag_arr, Sfab, TaggingParm::denerr, tagval);
          });
      }
      if (level < TaggingParm::max_dengrad_lev) {
        amrex::ParallelFor(
          tilebox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            tag_graderror(i, j, k, tag_arr, Sfab, TaggingParm::dengrad, tagval);
          });
      }

      // Tagging pressure
      S_derData.setVal<amrex::RunOn::Device>(0.0);
      erf_derpres(
        datbox, S_derData, ncp, Sfab.nComp(), S_data[mfi], geom, time, bc,
        level);
      if (level < TaggingParm::max_presserr_lev) {
        amrex::ParallelFor(
          tilebox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            tag_error(
              i, j, k, tag_arr, S_derarr, TaggingParm::presserr, tagval);
          });
      }
      if (level < TaggingParm::max_pressgrad_lev) {
        amrex::ParallelFor(
          tilebox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            tag_graderror(
              i, j, k, tag_arr, S_derarr, TaggingParm::pressgrad, tagval);
          });
      }

      // Tagging temperature
      S_derData.setVal<amrex::RunOn::Device>(0.0);
      erf_dertemp(
        datbox, S_derData, ncp, Sfab.nComp(), S_data[mfi], geom, time, bc,
        level);
      if (level < TaggingParm::max_temperr_lev) {
        amrex::ParallelFor(
          tilebox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            tag_error(i, j, k, tag_arr, S_derarr, TaggingParm::temperr, tagval);
          });
      }
      if (level < TaggingParm::max_tempgrad_lev) {
        amrex::ParallelFor(
          tilebox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            tag_graderror(
              i, j, k, tag_arr, S_derarr, TaggingParm::tempgrad, tagval);
          });
      }

      // Problem specific tagging
      const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx =
        geom.CellSizeArray();
      const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> prob_lo =
        geom.ProbLoArray();
      const auto captured_level = level;
      amrex::ParallelFor(
        tilebox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          set_problem_tags<ProblemTags>(
            i, j, k, tag_arr, Sfab, tagval, dx, prob_lo, time, captured_level);
        });

      // Now update the tags in the TagBox.
      // tag_arr.tags(itags, tilebox);
      // rho_eli.clear();
      // temp_eli.clear();
    }
  }
}

std::unique_ptr<amrex::MultiFab>
ERF::derive(const std::string& name, amrex::Real time, int ngrow)
{
  if (name == "x_velocity") 
  {
      MultiFab const& x_vel_on_face = get_new_data(X_Vel_Type);
      MultiFab const& y_vel_on_face = get_new_data(Y_Vel_Type);
      MultiFab const& z_vel_on_face = get_new_data(Z_Vel_Type);

      MultiFab ccvel(grids,dmap,3,0);

      average_face_to_cellcenter(ccvel,0,
          Array<const MultiFab*,3>{&x_vel_on_face, &y_vel_on_face, &z_vel_on_face});

      std::unique_ptr<MultiFab> derive_dat (new MultiFab(grids, dmap, 1, 0));
      MultiFab::Copy(*derive_dat, ccvel, 0, 0, 1, 0);

      return std::move(derive_dat);
  } 
  else if (name == "y_velocity") 
  {
      MultiFab const& x_vel_on_face = get_new_data(X_Vel_Type);
      MultiFab const& y_vel_on_face = get_new_data(Y_Vel_Type);
      MultiFab const& z_vel_on_face = get_new_data(Z_Vel_Type);

      MultiFab ccvel(grids,dmap,3,0);

      average_face_to_cellcenter(ccvel,0,
          Array<const MultiFab*,3>{&x_vel_on_face, &y_vel_on_face, &z_vel_on_face});

      std::unique_ptr<MultiFab> derive_dat (new MultiFab(grids, dmap, 1, 0));
      MultiFab::Copy(*derive_dat, ccvel, 1, 0, 1, 0);

      return std::move(derive_dat);
  } 
  else if (name == "z_velocity") 
  {
      MultiFab const& x_vel_on_face = get_new_data(X_Vel_Type);
      MultiFab const& y_vel_on_face = get_new_data(Y_Vel_Type);
      MultiFab const& z_vel_on_face = get_new_data(Z_Vel_Type);

      MultiFab ccvel(grids,dmap,3,0);

      average_face_to_cellcenter(ccvel,0,
          Array<const MultiFab*,3>{&x_vel_on_face, &y_vel_on_face, &z_vel_on_face});

      std::unique_ptr<MultiFab> derive_dat (new MultiFab(grids, dmap, 1, 0));
      MultiFab::Copy(*derive_dat, ccvel, 2, 0, 1, 0);

      return std::move(derive_dat);

  } else {
     return AmrLevel::derive(name, time, ngrow);
  } 
}

void
ERF::derive(
  const std::string& name, amrex::Real time, amrex::MultiFab& mf, int dcomp)
{
  {
    AmrLevel::derive(name, time, mf, dcomp);
  }
}

void
ERF::clear_prob()
{
  erf_prob_close();
}

void
ERF::init_transport()
{
  transport_init();
}

#ifdef ERF_USE_MASA
void
ERF::init_mms()
{
  if (!mms_initialized) {
    if (verbose && amrex::ParallelDescriptor::IOProcessor()) {
      amrex::Print() << "Initializing MMS" << std::endl;
    }
#ifdef ERF_USE_MASA
    masa_init("mms", masa_solution_name.c_str());
    masa_set_param("Cs", ERF::Cs);
    masa_set_param("CI", ERF::CI);
    masa_set_param("PrT", ERF::PrT);
#endif
    mms_initialized = true;
  }
}
#endif

amrex::Real
ERF::getCPUTime()
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

amrex::MultiFab&
ERF::build_fine_mask()
{
  // Mask for zeroing covered cells
  AMREX_ASSERT(level > 0);

  if (!fine_mask.empty())
    return fine_mask;

  const amrex::BoxArray& cba = parent->boxArray(level - 1);
  const amrex::DistributionMapping& cdm = parent->DistributionMap(level - 1);

  fine_mask.define(cba, cdm, 1, 0, amrex::MFInfo(), amrex::FArrayBoxFactory());
  fine_mask.setVal(1.0);

  amrex::BoxArray fba = parent->boxArray(level);
  amrex::iMultiFab ifine_mask = makeFineMask(cba, cdm, fba, crse_ratio, 1, 0);

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
  for (amrex::MFIter mfi(fine_mask, amrex::TilingIfNotGPU()); mfi.isValid();
       ++mfi) {
    auto& fab = fine_mask[mfi];
    auto& ifab = ifine_mask[mfi];
    const auto arr = fab.array();
    const auto iarr = ifab.array();
    amrex::ParallelFor(
      fab.box(), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
#ifdef _OPENMP
#pragma omp atomic write
#endif
        arr(i, j, k) = iarr(i, j, k);
      });
  }

  return fine_mask;
}

const amrex::iMultiFab*
ERF::build_interior_boundary_mask(int ng)
{
  for (int i = 0; i < ib_mask.size(); ++i) {
    if (ib_mask[i]->nGrow() == ng) {
      return ib_mask[i].get();
    }
  }

  //  If we got here, we need to build a new one
  if (ib_mask.size() == 0) {
    ib_mask.resize(0);
  }

  ib_mask.push_back(std::unique_ptr<amrex::iMultiFab>(new amrex::iMultiFab(
    grids, dmap, 1, ng, amrex::MFInfo(),
    amrex::DefaultFabFactory<amrex::IArrayBox>())));

  amrex::iMultiFab* imf = ib_mask.back().get();
  int ghost_covered_by_valid = 0;
  int other_cells =
    1; // uncovered ghost, valid, and outside domain cells are set to 1

  imf->BuildMask(
    geom.Domain(), geom.periodicity(), ghost_covered_by_valid, other_cells,
    other_cells, other_cells);

  return imf;
}
