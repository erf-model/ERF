#include <cstdio>

#include <AMReX_LevelBld.H>
#include <AMReX_ParmParse.H>
#include <AMReX_buildInfo.H>

#ifdef ERF_USE_MASA
#include <masa.h>
using namespace MASA;
#endif

#include "ERF.H"
#include "Derive.H"
#include "IndexDefines.H"
#include "prob.H"

//
// Components are:
//  Interior, Inflow, Outflow,  Symmetry,     SlipWall,     NoSlipWall, UserBC
//
static int scalar_bc[] = {INT_DIR,      EXT_DIR,      FOEXTRAP, REFLECT_EVEN,
                          REFLECT_EVEN, REFLECT_EVEN, EXT_DIR};

static int norm_vel_bc[] = {INT_DIR,     EXT_DIR,     FOEXTRAP, REFLECT_ODD,
                            REFLECT_ODD, REFLECT_ODD, EXT_DIR};

static int tang_vel_bc[] = {INT_DIR,      EXT_DIR,     FOEXTRAP, REFLECT_EVEN,
                            REFLECT_EVEN, REFLECT_ODD, EXT_DIR};

static void
set_scalar_bc(amrex::BCRec& bc, const amrex::BCRec& phys_bc)
{
  const int* lo_bc = phys_bc.lo();
  const int* hi_bc = phys_bc.hi();
  for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
    bc.setLo(dir, scalar_bc[lo_bc[dir]]);
    bc.setHi(dir, scalar_bc[hi_bc[dir]]);
  }
}

static void
set_x_vel_bc(amrex::BCRec& bc, const amrex::BCRec& phys_bc)
{
  const int* lo_bc = phys_bc.lo();
  const int* hi_bc = phys_bc.hi();
  AMREX_D_TERM(
    bc.setLo(0, norm_vel_bc[lo_bc[0]]); bc.setHi(0, norm_vel_bc[hi_bc[0]]);
    , bc.setLo(1, tang_vel_bc[lo_bc[1]]); bc.setHi(1, tang_vel_bc[hi_bc[1]]);
    , bc.setLo(2, tang_vel_bc[lo_bc[2]]); bc.setHi(2, tang_vel_bc[hi_bc[2]]););
}

static void
set_y_vel_bc(amrex::BCRec& bc, const amrex::BCRec& phys_bc)
{
  const int* lo_bc = phys_bc.lo();
  const int* hi_bc = phys_bc.hi();
  AMREX_D_TERM(
    bc.setLo(0, tang_vel_bc[lo_bc[0]]); bc.setHi(0, tang_vel_bc[hi_bc[0]]);
    , bc.setLo(1, norm_vel_bc[lo_bc[1]]); bc.setHi(1, norm_vel_bc[hi_bc[1]]);
    , bc.setLo(2, tang_vel_bc[lo_bc[2]]); bc.setHi(2, tang_vel_bc[hi_bc[2]]););
}

static void
set_z_vel_bc(amrex::BCRec& bc, const amrex::BCRec& phys_bc)
{
  const int* lo_bc = phys_bc.lo();
  const int* hi_bc = phys_bc.hi();
  AMREX_D_TERM(
    bc.setLo(0, tang_vel_bc[lo_bc[0]]); bc.setHi(0, tang_vel_bc[hi_bc[0]]);
    , bc.setLo(1, tang_vel_bc[lo_bc[1]]); bc.setHi(1, tang_vel_bc[hi_bc[1]]);
    , bc.setLo(2, norm_vel_bc[lo_bc[2]]); bc.setHi(2, norm_vel_bc[hi_bc[2]]););
}

void
ERF::variableSetUp()
{
  // ERF::variableSetUp is called in the constructor of Amr.cpp, so
  // it should get called every time we start or restart a job

  // initialize the start time for our CPU-time tracker
  startCPUTime = amrex::ParallelDescriptor::second();

  // Output the git commit hashes used to build the executable.

  if (amrex::ParallelDescriptor::IOProcessor()) {
    const char* erf_hash = amrex::buildInfoGetGitHash(1);
    const char* amrex_hash = amrex::buildInfoGetGitHash(2);
    const char* buildgithash = amrex::buildInfoGetBuildGitHash();
    const char* buildgitname = amrex::buildInfoGetBuildGitName();

    if (strlen(erf_hash) > 0) {
      amrex::Print() << "\n"
                     << "ERF git hash: " << erf_hash << "\n";
    }
    if (strlen(amrex_hash) > 0) {
      amrex::Print() << "AMReX git hash: " << amrex_hash << "\n";
    }
    if (strlen(buildgithash) > 0) {
      amrex::Print() << buildgitname << " git hash: " << buildgithash << "\n";
    }

    amrex::Print() << "\n";
  }

  AMREX_ASSERT(desc_lst.size() == 0);

  // Get options, set phys_bc
  read_params();

  init_transport();

  indxmap::init();

#ifdef ERF_USE_MASA
  if (do_mms) {
    init_mms();
  }
#endif

  //
  // Set number of state variables and pointers to components
  //

#ifdef NUM_ADV
  NumAdv = NUM_ADV;
#else
  NumAdv = 0;
#endif

  if (NumAdv > 0) {
    FirstAdv = Theta_comp+1;
  }

  int dm = AMREX_SPACEDIM;

  amrex::Vector<amrex::Real> center(AMREX_SPACEDIM, 0.0);
  amrex::ParmParse ppc("erf");
  ppc.queryarr("center", center, 0, AMREX_SPACEDIM);

  amrex::Interpolater* interp;

  if (state_interp_order == 0) {
    interp = &amrex::pc_interp;
  } else {
    if (lin_limit_state_interp == 1) {
      interp = &amrex::lincc_interp;
    } else {
      interp = &amrex::cell_cons_interp;
    }
  }

  // Note that the default is state_data_extrap = false,
  // store_in_checkpoint = true.  We only need to put these in
  // explicitly if we want to do something different,
  // like not store the state data in a checkpoint directory
  bool state_data_extrap = false;
  bool store_in_checkpoint;

  int ngrow_state = 2;

  // NVAR is currently set to 2 in IndexDefines.H
  store_in_checkpoint = true;
  desc_lst.addDescriptor(
    State_Type, amrex::IndexType::TheCellType(), amrex::StateDescriptor::Point,
    ngrow_state, NVAR, interp, state_data_extrap, store_in_checkpoint);

  amrex::Vector<amrex::BCRec> bcs(NVAR);
  amrex::Vector<std::string> name(NVAR);

  amrex::BCRec bc;
  int cnt = 0;
  set_scalar_bc(bc, phys_bc);
  bcs[cnt] = bc;
  name[cnt] = "density";
  cnt++;
  set_scalar_bc(bc, phys_bc);
  bcs[cnt] = bc;
  name[cnt] = "theta";

  for (int i = 0; i < NumAdv; ++i) {
    char buf[64];
    sprintf(buf, "adv_%d", i);
    cnt++;
    set_scalar_bc(bc, phys_bc);
    bcs[cnt] = bc;
    name[cnt] = std::string(buf);
  }

  amrex::StateDescriptor::BndryFunc bndryfunc1(erf_bcfill_hyp);
  bndryfunc1.setRunOnGPU(true);

  desc_lst.setComponent(State_Type, Density_comp, name, bcs, bndryfunc1);

  if (do_mol_load_balance) {
    desc_lst.addDescriptor(
      Work_Estimate_Type, amrex::IndexType::TheCellType(),
      amrex::StateDescriptor::Point, 0, 1, &amrex::pc_interp);
    // Because we use piecewise constant interpolation, we do not use bc and
    // BndryFunc.
    desc_lst.setComponent(
      Work_Estimate_Type, 0, "WorkEstimate", bc,
      amrex::StateDescriptor::BndryFunc(erf_nullfill));
  }

  // 
  // Create face-based StateData for velocities on each face,
  //     and for all conserved quantities as well
  //     (unnecessary right now but convenient)
  // 
  store_in_checkpoint = true;
  amrex::IndexType xface(amrex::IntVect(1,0,0));
  desc_lst.addDescriptor(X_Vel_Type, xface,
                         amrex::StateDescriptor::Point, 1, 1,
                         interp, state_data_extrap,
                         store_in_checkpoint);

  amrex::IndexType yface(amrex::IntVect(0,1,0));
  desc_lst.addDescriptor(Y_Vel_Type, yface,
                         amrex::StateDescriptor::Point, 1, 1,
                         interp, state_data_extrap,
                         store_in_checkpoint);

  amrex::IndexType zface(amrex::IntVect(0,0,1));
  desc_lst.addDescriptor(Z_Vel_Type, zface,
                         amrex::StateDescriptor::Point, 1, 1,
                         interp, state_data_extrap,
                         store_in_checkpoint);

  num_state_type = desc_lst.size();

  //
  // DEFINE DERIVED QUANTITIES
  //
  // Pressure
  //
  derive_lst.add(
    "pressure", amrex::IndexType::TheCellType(), 1, pc_derpres, the_same_box);
  derive_lst.addComponent("pressure", desc_lst, State_Type, Density_comp, NVAR);

  //
  // Temperature
  //
  derive_lst.add(
    "temp", amrex::IndexType::TheCellType(), 1, pc_dertemp, the_same_box);
  derive_lst.addComponent("temp", desc_lst, State_Type, Density_comp, NVAR);

  //
  // Velocities
  //

  // This calculcates cell-centered x-velocity from x-face-centered values
  // Note that we don't need to use "addComponent" because we are are computing
  //      these in ERF itself, not calling the AMReX derive routines
  derive_lst.add(
    "x_velocity", amrex::IndexType::TheCellType(), 1, pc_dernull, the_same_box);

  // This calculcates cell-centered y-velocity from y-face-centered values
  derive_lst.add(
    "y_velocity", amrex::IndexType::TheCellType(), 1, pc_dernull, the_same_box);

  // This calculcates cell-centered z-velocity from z-face-centered values
  derive_lst.add(
    "z_velocity", amrex::IndexType::TheCellType(), 1, pc_dernull, the_same_box);

  // Problem-specific derives
  add_problem_derives<ProblemDerives>(derive_lst, desc_lst);

  // Set list of active sources
  set_active_sources();
}

void
ERF::set_active_sources()
{
  // optional external source
  if (add_ext_src == 1) {
    src_list.push_back(ext_src);
  }

  // optional forcing source
  if (add_forcing_src == 1) {
    src_list.push_back(forcing_src);
  }
}
