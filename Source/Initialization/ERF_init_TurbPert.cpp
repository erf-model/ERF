// May 28th, 2024
// DUSTIN MA

#include <ERF.H>
#include <ERF_Constants.H>
#include <TileNoZ.H>
#include <prob_common.H>

using namespace amrex;

void
ERF::init_TurbPert_updateTime (int lev, TurbulentPerturbation& turbPert)
{
    auto& lev_new = vars_new[lev];

    // Grabing data from velocity field
    MultiFab xvel_data(lev_new[Vars::xvel].boxArray(), lev_new[Vars::xvel].DistributionMap(), 1, lev_new[Vars::xvel].nGrowVect());
    MultiFab yvel_data(lev_new[Vars::yvel].boxArray(), lev_new[Vars::yvel].DistributionMap(), 1, lev_new[Vars::yvel].nGrowVect());
    MultiFab::Copy (xvel_data, lev_new[Vars::xvel], 0, 0, 1, lev_new[Vars::xvel].nGrowVect());
    MultiFab::Copy (yvel_data, lev_new[Vars::yvel], 0, 0, 1, lev_new[Vars::yvel].nGrowVect());

    for (MFIter mfi(lev_new[Vars::cons], TileNoZ()); mfi.isValid(); ++mfi) {
        const Box &bx  = mfi.tilebox();
        const Box &xbx = mfi.tilebox(IntVect(1,0,0));
        const Box &ybx = mfi.tilebox(IntVect(0,1,0));

        const auto &xvel_data_arr = xvel_data.array(mfi);
        const auto &yvel_data_arr = yvel_data.array(mfi);

        // Computing initial perturbation frequency
        turbPert.calc_tpi_updateTime(lev, 0., bx, xvel_data_arr, yvel_data_arr);
    }
    Print() << "Turbulent perturbation update time initialized\n";
}

// Calculate the perturbation region amplitude.
// This function heavily emmulates the ERF::init_custom ()
void
ERF::init_TurbPert_amplitude (int lev, TurbulentPerturbation& turbPert)
{
    auto& lev_new = vars_new[lev];

    MultiFab r_hse(base_state[lev], make_alias, 0, 1); // r_0 is first  component
    MultiFab p_hse(base_state[lev], make_alias, 1, 1); // p_0 is second component

    MultiFab cons_pert(lev_new[Vars::cons].boxArray(), lev_new[Vars::cons].DistributionMap(),
                       lev_new[Vars::cons].nComp()   , lev_new[Vars::cons].nGrow());
    MultiFab xvel_pert(lev_new[Vars::xvel].boxArray(), lev_new[Vars::xvel].DistributionMap(), 1, lev_new[Vars::xvel].nGrowVect());
    MultiFab yvel_pert(lev_new[Vars::yvel].boxArray(), lev_new[Vars::yvel].DistributionMap(), 1, lev_new[Vars::yvel].nGrowVect());
    MultiFab zvel_pert(lev_new[Vars::zvel].boxArray(), lev_new[Vars::zvel].DistributionMap(), 1, lev_new[Vars::zvel].nGrowVect());

    // Default perturbations to zero
    cons_pert.setVal(0.);
    xvel_pert.setVal(0.);
    yvel_pert.setVal(0.);
    zvel_pert.setVal(0.);

    int fix_random_seed = 0;
    ParmParse pp("erf"); pp.query("fix_random_seed", fix_random_seed);
    // Note that the value of 1024UL is not significant -- the point here is just to set the
    //     same seed for all MPI processes for the purpose of regression testing
    if (fix_random_seed) {
        Print() << "Fixing the random seed" << std::endl;
        InitRandom(1024UL);
    }

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(lev_new[Vars::cons], TileNoZ()); mfi.isValid(); ++mfi) {
        const Box &bx  = mfi.tilebox();
        const Box &xbx = mfi.tilebox(IntVect(1,0,0));
        const Box &ybx = mfi.tilebox(IntVect(0,1,0));
        const Box &zbx = mfi.tilebox(IntVect(0,0,1));

        // Perturbation on to different components
        const auto &cons_pert_arr = cons_pert.array(mfi);
        const auto &xvel_pert_arr = xvel_pert.array(mfi);
        const auto &yvel_pert_arr = yvel_pert.array(mfi);
        const auto &zvel_pert_arr = zvel_pert.array(mfi);

        Array4<Real const> cons_arr = lev_new[Vars::cons].const_array(mfi);
        Array4<Real const> z_nd_arr = (solverChoice.use_terrain) ? z_phys_nd[lev]->const_array(mfi) : Array4<Real const>{};
        Array4<Real const> z_cc_arr = (solverChoice.use_terrain) ? z_phys_cc[lev]->const_array(mfi) : Array4<Real const>{};

        Array4<Real const> mf_m     = mapfac_m[lev]->array(mfi);
        Array4<Real const> mf_u     = mapfac_m[lev]->array(mfi);
        Array4<Real const> mf_v     = mapfac_m[lev]->array(mfi);

        Array4<Real> r_hse_arr = r_hse.array(mfi);
        Array4<Real> p_hse_arr = p_hse.array(mfi);

        // Computing perturbation amplitude
        turbPert.apply_tpi(lev,  bx, RhoTheta_comp, cons_pert_arr); // Using boxArray

        prob->init_custom_pert(bx, xbx, ybx, zbx, cons_arr, cons_pert_arr,
                               xvel_pert_arr, yvel_pert_arr, zvel_pert_arr,
                               r_hse_arr, p_hse_arr, z_nd_arr, z_cc_arr,
                               geom[lev].data(), mf_m, mf_u, mf_v,
                               solverChoice);
    } // mfi
    MultiFab::Add(lev_new[Vars::cons], cons_pert, RhoTheta_comp, RhoTheta_comp, 1, cons_pert.nGrow());

    Print() << "Perturbation region amplitude initialized using : "
    #ifdef RANDOM_PERTURB
            << "Random number perturbation amplitude" 
    #else
            << "Artificial number index fill"
    #endif
            << "\n";
}
