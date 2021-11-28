/**
 * \file ERF.cpp
 */

#ifdef _OPENMP
#include <omp.h>
#endif

#include <AMReX_Vector.H>
#include <AMReX_Utility.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_VisMF.H>
#include <AMReX_TagBox.H>
#include <AMReX_ParmParse.H>

#include "ERF.H"
#include "ERF_Constants.H"
#include "Derive.H"
#include "EOS.H"
#include "Tagging.H"
#include "IndexDefines.H"

using namespace amrex;

bool ERF::signalStopJob = false;
bool ERF::dump_old = false;
int ERF::verbose = 0;

SolverChoice ERF::solverChoice;
ABLFieldInit ERF::ablinit;

amrex::Real ERF::cfl         = 0.8;
amrex::Real ERF::init_shrink = 1.0;
amrex::Real ERF::change_max  = 1.1;
amrex::Real ERF::initial_dt  = -1.0;
amrex::Real ERF::fixed_dt    = -1.0;
amrex::Real ERF::max_dt      =  1.e20;
amrex::Real ERF::dt_cutoff   =  0.0;

std::string ERF::coupling_type = "OneWay";
int         ERF::do_reflux     = 0;
int         ERF::do_avg_down   = 0;
int         ERF::sum_interval  = -1;
amrex::Real ERF::sum_per       = -1.0;

bool        ERF::init_abl      = false;

// Create dens_hse and pres_hse with one ghost cell
int ERF::ng_dens_hse = 1;
int ERF::ng_pres_hse = 1;

amrex::Vector<std::unique_ptr<phys_bcs::BCBase> > ERF::bc_recs(AMREX_SPACEDIM*2);
int ERF::NumAdv = 0;
int ERF::FirstAdv = -1;

bool ERF::lo_z_is_no_slip = false;
bool ERF::hi_z_is_no_slip = false;

amrex::Vector<int> ERF::src_list;

amrex::Vector<amrex::Vector<amrex::Real> > ERF::h_dens_hse(0);
amrex::Vector<amrex::Vector<amrex::Real> > ERF::h_pres_hse(0);
amrex::Vector<amrex::Gpu::DeviceVector<amrex::Real> > ERF::d_dens_hse(0);
amrex::Vector<amrex::Gpu::DeviceVector<amrex::Real> > ERF::d_pres_hse(0);

amrex::Vector<amrex::Vector<amrex::Real> > ERF::h_rayleigh_tau(0);
amrex::Vector<amrex::Vector<amrex::Real> > ERF::h_rayleigh_ubar(0);
amrex::Vector<amrex::Vector<amrex::Real> > ERF::h_rayleigh_vbar(0);
amrex::Vector<amrex::Vector<amrex::Real> > ERF::h_rayleigh_thetabar(0);
amrex::Vector<amrex::Gpu::DeviceVector<amrex::Real> > ERF::d_rayleigh_tau(0);
amrex::Vector<amrex::Gpu::DeviceVector<amrex::Real> > ERF::d_rayleigh_ubar(0);
amrex::Vector<amrex::Gpu::DeviceVector<amrex::Real> > ERF::d_rayleigh_vbar(0);
amrex::Vector<amrex::Gpu::DeviceVector<amrex::Real> > ERF::d_rayleigh_thetabar(0);

// this will be reset upon restart
amrex::Real ERF::previousCPUTimeUsed = 0.0;
amrex::Real ERF::startCPUTime = 0.0;
int ERF::num_state_type = 0;

amrex::Vector<std::string> BCNames = {"xlo", "ylo", "zlo", "xhi", "yhi", "zhi"};

template<int IDIR, math_bcs::BCBound Bound> std::string getBCName() {
    return BCNames[IDIR + 3*Bound];
}

template<int DIM, math_bcs::BCBound Bound>
std::unique_ptr<phys_bcs::BCBase>
ERF::initialize_bcs(const std::string& bc_char) {
  if (!bc_char.compare("Interior")) {
    std::unique_ptr<phys_bcs::BCBase> bc_rec(new phys_bcs::BCInterior());
    return std::move(bc_rec);
  } else if (!bc_char.compare("Hard")) {
    std::unique_ptr<phys_bcs::BCBase> bc_rec(new phys_bcs::BCDummy());
    return std::move(bc_rec);
  } else if (!bc_char.compare("FOExtrap")) {
    std::unique_ptr<phys_bcs::BCBase> bc_rec(new phys_bcs::BCDummy());
    return std::move(bc_rec);
  } else if (!bc_char.compare("Symmetry")) {
    std::unique_ptr<phys_bcs::BCBase> bc_rec(new phys_bcs::BCDummy());
    return std::move(bc_rec);
  } else if (!bc_char.compare("Dirichlet")) {
    amrex::ParmParse bcinp(getBCName<DIM, Bound>());
    amrex::Vector<amrex::Real> uvec;
    bcinp.getarr("velocity", uvec, 0, AMREX_SPACEDIM);
    amrex::Print() << "Dirichlet selected for DIM=" << DIM
        << " lower/upper=" << Bound
        << " fixedvalue=" << uvec[0] << " " << uvec[1] << " " << uvec[2]
        << std::endl;
    std::unique_ptr<phys_bcs::BCBase> bc_rec(new phys_bcs::BCDirichlet<DIM, Bound>(uvec));
    return std::move(bc_rec);
  } else if (!bc_char.compare("SlipWall")) {
    std::unique_ptr<phys_bcs::BCBase> bc_rec(new phys_bcs::BCSlipWall<DIM, Bound>());
    return std::move(bc_rec);
  } else if (!bc_char.compare("NoSlipWall")) {
    std::unique_ptr<phys_bcs::BCBase> bc_rec(new phys_bcs::BCNoSlipWall<DIM, Bound>());
    return std::move(bc_rec);
  } else if (!bc_char.compare("SimSlipWall")) {
    std::unique_ptr<phys_bcs::BCBase> bc_rec(new phys_bcs::BCSimSlipWall<DIM, Bound>());
    return std::move(bc_rec);
  } else {
    amrex::Abort("Wrong boundary condition word, please use: "
                 "Interior, Dirichlet, SimSlipWall, Symmetry, SlipWall, NoSlipWall");
    return NULL;
  }
}

void
ERF::read_params()
{
  static bool read_params_done = false;

  if (read_params_done)
    return;

  read_params_done = true;

  amrex::ParmParse pp("erf");

  pp.query("v", verbose);
  pp.query("sum_interval", sum_interval);
  pp.query("dump_old", dump_old);

  pp.query("coupling_type",coupling_type);
  if (coupling_type == "OneWay")
  {
      do_reflux = 0;
      do_avg_down = 0;
  }
  else if (coupling_type == "TwoWay")
  {
      do_reflux = 1;
      do_avg_down = 1;
  } else {
      amrex::Error("Unknown coupling type");
  }

  // This defaults to false; is only true if we want to call ABLFieldInit::init_params
  pp.query("init_abl", init_abl);

  // Time step controls
  pp.query("cfl", cfl);
  pp.query("init_shrink", init_shrink);
  pp.query("change_max", change_max);
  pp.query("initial_dt", initial_dt);
  pp.query("fixed_dt", fixed_dt);
  pp.query("max_dt", max_dt);
  pp.query("dt_cutoff", dt_cutoff);

  // Get boundary conditions
  amrex::Vector<std::string> lo_bc_char(AMREX_SPACEDIM);
  amrex::Vector<std::string> hi_bc_char(AMREX_SPACEDIM);
  pp.getarr("lo_bc", lo_bc_char, 0, AMREX_SPACEDIM);
  pp.getarr("hi_bc", hi_bc_char, 0, AMREX_SPACEDIM);

  if (lo_bc_char[0] == "NoSlipWall" || hi_bc_char[0] == "NoSlipWall" ||
      lo_bc_char[1] == "NoSlipWall" || hi_bc_char[1] == "NoSlipWall")
      amrex::Error("No-slip wall only allowed on z-faces");

  bc_recs[0] = ERF::initialize_bcs<0, math_bcs::BCBound::lower>(lo_bc_char[0]);
  bc_recs[1] = ERF::initialize_bcs<0, math_bcs::BCBound::upper>(hi_bc_char[0]);
  bc_recs[2] = ERF::initialize_bcs<1, math_bcs::BCBound::lower>(lo_bc_char[1]);
  bc_recs[3] = ERF::initialize_bcs<1, math_bcs::BCBound::upper>(hi_bc_char[1]);
  bc_recs[4] = ERF::initialize_bcs<2, math_bcs::BCBound::lower>(lo_bc_char[2]);
  bc_recs[5] = ERF::initialize_bcs<2, math_bcs::BCBound::upper>(hi_bc_char[2]);

  if (lo_bc_char[2] == "NoSlipWall") lo_z_is_no_slip = true;
  if (hi_bc_char[2] == "NoSlipWall") hi_z_is_no_slip = true;

  //
  // Check bc_recs against possible periodic geometry
  // if periodic, must have internal BC marked.
  //
  //
  // Do idiot check.  Periodic means interior in those directions.
  //
  for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
    if (amrex::DefaultGeometry().isPeriodic(dir)) {
      if (
        !(bc_recs[2*dir]->isInterior()) && amrex::ParallelDescriptor::IOProcessor()) {
        std::cerr << "ERF::read_params:periodic in direction " << dir
                  << " but low BC is not Interior\n";
        amrex::Error();
      }
      if (
        !(bc_recs[2*dir+1]->isInterior()) && amrex::ParallelDescriptor::IOProcessor()) {
        std::cerr << "ERF::read_params:periodic in direction " << dir
                  << " but high BC is not Interior\n";
        amrex::Error();
      }
    } else {
      //
      // Do idiot check. If not periodic, should not be interior.
      //
      if (bc_recs[2*dir]->isInterior() && amrex::ParallelDescriptor::IOProcessor()) {
        std::cerr << "ERF::read_params:interior bc in direction " << dir
                  << " but not periodic\n";
        amrex::Error();
      }
      if (bc_recs[2*dir+1]->isInterior() && amrex::ParallelDescriptor::IOProcessor()) {
        std::cerr << "ERF::read_params:interior bc in direction " << dir
                  << " but not periodic\n";
        amrex::Error();
      }
    }
  }

  // Sanity check
  if (cfl <= 0.0 || cfl > 1.0) {
    amrex::Error("Invalid CFL factor; must be between zero and one.");
  }

  if (max_dt < fixed_dt) {
    amrex::Error("Cannot have max_dt < fixed_dt");
  }

  // Read tagging parameters
  read_tagging_params();

  solverChoice.init_params();
}

ERF::ERF()
  : io_mgr(new IOManager(*this))
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
    io_mgr(new IOManager(*this))
{
  buildMetrics();

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
//  const int ngrd = grids.size();

//  const amrex::Real* dx = geom.CellSize();

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
}

void
ERF::setTimeLevel(amrex::Real time, amrex::Real dt_old, amrex::Real dt_new)
{
  AmrLevel::setTimeLevel(time, dt_old, dt_new);
}

void
ERF::initData()
{
  BL_PROFILE("ERF::initData()");

  amrex::MultiFab& S_new = get_new_data(State_Type);
  amrex::MultiFab& U_new = get_new_data(X_Vel_Type);
  amrex::MultiFab& V_new = get_new_data(Y_Vel_Type);
  amrex::MultiFab& W_new = get_new_data(Z_Vel_Type);

  if (level == 0) {

    init1DArrays();

    if (init_abl)
    {
        ablinit.init_params();
    }
  }

  initDataProb(S_new, U_new, V_new, W_new);
}

void
ERF::init1DArrays()
{
    AMREX_ALWAYS_ASSERT(level == 0);
    //
    // Setup Base State Arrays
    //
    const int max_level = parent->maxLevel();
    h_dens_hse.resize(max_level+1, amrex::Vector<Real>(0));
    h_pres_hse.resize(max_level+1, amrex::Vector<Real>(0));
    d_dens_hse.resize(max_level+1, amrex::Gpu::DeviceVector<Real>(0));
    d_pres_hse.resize(max_level+1, amrex::Gpu::DeviceVector<Real>(0));

    if (solverChoice.use_rayleigh_damping)
    {
        h_rayleigh_tau.resize(max_level+1, amrex::Vector<Real>(0));
        h_rayleigh_ubar.resize(max_level+1, amrex::Vector<Real>(0));
        h_rayleigh_vbar.resize(max_level+1, amrex::Vector<Real>(0));
        h_rayleigh_thetabar.resize(max_level+1, amrex::Vector<Real>(0));
        d_rayleigh_tau.resize(max_level+1, amrex::Gpu::DeviceVector<Real>(0));
        d_rayleigh_ubar.resize(max_level+1, amrex::Gpu::DeviceVector<Real>(0));
        d_rayleigh_vbar.resize(max_level+1, amrex::Gpu::DeviceVector<Real>(0));
        d_rayleigh_thetabar.resize(max_level+1, amrex::Gpu::DeviceVector<Real>(0));
    }

    // Initialize the hydrostatic and Rayleigh damping arrays at all levels
    //     and set from problem-dependent input
    for (int lev(0); lev <= max_level; ++lev) {

      const int zlen_dens = parent->Geom(lev).Domain().length(2) + 2*ng_dens_hse;
      h_dens_hse[lev].resize(zlen_dens, 0.0_rt);
      d_dens_hse[lev].resize(zlen_dens, 0.0_rt);

      const int zlen_pres = parent->Geom(lev).Domain().length(2) + 2*ng_pres_hse;
      h_pres_hse[lev].resize(zlen_pres, p_0);
      d_pres_hse[lev].resize(zlen_pres, p_0);

      getLevel(lev).initHSE();

      if (solverChoice.use_rayleigh_damping)
      {
          const int zlen_rayleigh = parent->Geom(lev).Domain().length(2);
          h_rayleigh_tau[lev].resize(zlen_rayleigh, 0.0_rt);
          d_rayleigh_tau[lev].resize(zlen_rayleigh, 0.0_rt);
          h_rayleigh_ubar[lev].resize(zlen_rayleigh, 0.0_rt);
          d_rayleigh_ubar[lev].resize(zlen_rayleigh, 0.0_rt);
          h_rayleigh_vbar[lev].resize(zlen_rayleigh, 0.0_rt);
          d_rayleigh_vbar[lev].resize(zlen_rayleigh, 0.0_rt);
          h_rayleigh_thetabar[lev].resize(zlen_rayleigh, 0.0_rt);
          d_rayleigh_thetabar[lev].resize(zlen_rayleigh, 0.0_rt);

          getLevel(lev).initRayleigh();
      }
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
  amrex::MultiFab& U_new = get_new_data(X_Vel_Type);
  amrex::MultiFab& V_new = get_new_data(Y_Vel_Type);
  amrex::MultiFab& W_new = get_new_data(Z_Vel_Type);

  FillPatch(old, S_new, 0, cur_time, State_Type, 0, NVAR);
  FillPatch(old, U_new, 0, cur_time, X_Vel_Type, 0, 1);
  FillPatch(old, V_new, 0, cur_time, Y_Vel_Type, 0, 1);
  FillPatch(old, W_new, 0, cur_time, Z_Vel_Type, 0, 1);
}

/**
 * This version inits the data on a new level that did not
 * exist before regridding.
 */
void
ERF::init()
{
  BL_PROFILE("ERF::init()");

  amrex::Real dt = parent->dtLevel(level);
  amrex::Real cur_time = getLevel(level - 1).state[State_Type].curTime();
  amrex::Real prev_time = getLevel(level - 1).state[State_Type].prevTime();

  amrex::Real dt_old =
    (cur_time - prev_time) / (amrex::Real)parent->MaxRefRatio(level - 1);

  setTimeLevel(cur_time, dt_old, dt);
  amrex::MultiFab& S_new = get_new_data(State_Type);
  FillCoarsePatch(S_new, 0, cur_time, State_Type, 0, NVAR);
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
               const amrex::Real rho      = s(i, j, k, Rho_comp);
               const amrex::Real rhotheta = s(i, j, k, RhoTheta_comp);

               amrex::Real pressure = getPgivenRTh(rhotheta);
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
//  const int ncycle = parent->nCycle(level);

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
//    amrex::MultiFab& S_new_crse = get_new_data(State_Type);
    // clean_state(S_new_crse);
  }

  // Re-compute temperature after all the other updates.
//  amrex::MultiFab& S_new = get_new_data(State_Type);
//  int ng_pts = 0;

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

  if (level == 0) {
    init1DArrays();
  }

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
  const char tagval = amrex::TagBox::SET;

  amrex::MultiFab S_data(
    get_new_data(State_Type).boxArray(),
    get_new_data(State_Type).DistributionMap(), NVAR, 1);

  // Static refinement of a specified region
  if (tparm.tag_region)
  {
      Real xlo = (tparm.region_lo[0] - geom.ProbLo()[0]);
      Real ylo = (tparm.region_lo[1] - geom.ProbLo()[0]);
      Real xhi = (tparm.region_hi[0] - geom.ProbLo()[1]);
      Real yhi = (tparm.region_hi[1] - geom.ProbLo()[1]);
      Real zlo = (tparm.region_lo[2] - geom.ProbLo()[2]);
      Real zhi = (tparm.region_hi[2] - geom.ProbLo()[2]);

      const Real l_dx = geom.CellSize(0);
      const Real l_dy = geom.CellSize(1);
      const Real l_dz = geom.CellSize(2);

      for (amrex::MFIter mfi(S_data, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
         const amrex::Box& bx = mfi.tilebox();
         auto tag_arr = tags.array(mfi);

         amrex::ParallelFor(bx,
         [xlo, xhi, ylo, yhi, zlo, zhi, l_dx, l_dy, l_dz,tagval, tag_arr]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
              Real x = (i+0.5)*l_dx;
              Real y = (j+0.5)*l_dy;
              Real z = (k+0.5)*l_dz;

              // Tag if we are inside the specified box
              if (x >= xlo && x <= xhi && y >= ylo && y <= yhi && z >= zlo && z <= zhi)
              {
                 tag_arr(i,j,k) = tagval;
              }
         });
      }

  } else {

  const amrex::Real cur_time = state[State_Type].curTime();
  FillPatch(
    *this, S_data, S_data.nGrow(), cur_time, State_Type, Rho_comp, NVAR, 0);

  amrex::Vector<amrex::BCRec> bcs(NVAR);
//  const char clearval = amrex::TagBox::CLEAR;

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
      if (level < tparm.max_denerr_lev) {
        amrex::ParallelFor(
          tilebox, [=, denerr=tparm.denerr] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            tag_error(i, j, k, tag_arr, Sfab, denerr, tagval);
          });
      }
      if (level < tparm.max_dengrad_lev) {
        amrex::ParallelFor(
          tilebox, [=, dengrad=tparm.dengrad] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            tag_graderror(i, j, k, tag_arr, Sfab, dengrad, tagval);
          });
      }

      // Tagging pressure
      S_derData.setVal<amrex::RunOn::Device>(0.0);
      erf_derpres(
        datbox, S_derData, ncp, Sfab.nComp(), S_data[mfi], geom, time, bc,
        level);
      if (level < tparm.max_presserr_lev) {
        amrex::ParallelFor(
          tilebox, [=, presserr=tparm.presserr] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            tag_error(
              i, j, k, tag_arr, S_derarr, presserr, tagval);
          });
      }
      if (level < tparm.max_pressgrad_lev) {
        amrex::ParallelFor(
          tilebox, [=, pressgrad=tparm.pressgrad] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            tag_graderror(
              i, j, k, tag_arr, S_derarr, pressgrad, tagval);
          });
      }

      // Tagging temperature
      S_derData.setVal<amrex::RunOn::Device>(0.0);
      erf_dertemp(
        datbox, S_derData, ncp, Sfab.nComp(), S_data[mfi], geom, time, bc,
        level);
      if (level < tparm.max_temperr_lev) {
        amrex::ParallelFor(
          tilebox, [=, temperr=tparm.temperr] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            tag_error(i, j, k, tag_arr, S_derarr, temperr, tagval);
          });
      }
      if (level < tparm.max_tempgrad_lev) {
        amrex::ParallelFor(
          tilebox, [=, tempgrad=tparm.tempgrad] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            tag_graderror(
              i, j, k, tag_arr, S_derarr, tempgrad, tagval);
          });
      }

      // Now update the tags in the TagBox.
      // tag_arr.tags(itags, tilebox);
    }
  }
  } // end of !tag_region
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

void ERF::restart(amrex::Amr& papa, istream& is, bool bReadSpecial)
{
#ifdef ERF_USE_NETCDF
  io_mgr->ncrestart(papa, is, bReadSpecial);
#else
  AmrLevel::restart(papa, is, bReadSpecial);
  io_mgr->restart(papa, is, bReadSpecial);
#endif
}

void ERF::set_state_in_checkpoint(amrex::Vector<int>& state_in_checkpoint)
{
  io_mgr->set_state_in_checkpoint(state_in_checkpoint);
}

void ERF::checkPoint(const std::string& dir, std::ostream& os,
                     amrex::VisMF::How how, bool dump_old_default)
{
#ifdef ERF_USE_NETCDF
  io_mgr->NCWriteCheckpointFile (dir, os, dump_old);
#else
  AmrLevel::checkPoint(dir, os, how, dump_old);
  io_mgr->checkPoint(dir, os, how, dump_old_default);
#endif
}

void ERF::setPlotVariables()
{
  amrex::AmrLevel::setPlotVariables();
  io_mgr->setPlotVariables();
}

void ERF::writeJobInfo(const std::string& dir)
{
  io_mgr->writeJobInfo(dir);
}

void ERF::writePlotFile(const std::string& dir, ostream& os, amrex::VisMF::How how)
{
#ifdef ERF_USE_NETCDF
  io_mgr->writeNCPlotFile(dir, os);
#else
  io_mgr->writePlotFile(dir, os, how);
#endif
}

void ERF::writeSmallPlotFile(const std::string& dir, ostream& os, amrex::VisMF::How how)
{
  io_mgr->writeSmallPlotFile(dir, os, how);
}

