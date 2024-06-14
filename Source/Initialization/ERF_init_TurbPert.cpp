// May 28th, 2024
// DUSTIN MA

#include <ERF.H>
#include <ERF_Constants.H>
#include <TileNoZ.H>
#include <prob_common.H>

using namespace amrex;

void
ERF::TurbPert_constants(const int lev)
{
    prob->init_turbPert_const(turbPert);
}

#define USE_SLAB_AVERAGE
//#define USE_VOLUME_AVERAGE

void
ERF::TurbPert_update (const int lev, const Real dt, TurbulentPerturbation& turbPert)
{
    // Grabing data from velocity field
    auto& lev_new = vars_new[lev];

    MultiFab xvel_data(lev_new[Vars::xvel].boxArray(), lev_new[Vars::xvel].DistributionMap(), 1, lev_new[Vars::xvel].nGrowVect());
    MultiFab yvel_data(lev_new[Vars::yvel].boxArray(), lev_new[Vars::yvel].DistributionMap(), 1, lev_new[Vars::yvel].nGrowVect());
    MultiFab::Copy (xvel_data, lev_new[Vars::xvel], 0, 0, 1, lev_new[Vars::xvel].nGrowVect());
    MultiFab::Copy (yvel_data, lev_new[Vars::yvel], 0, 0, 1, lev_new[Vars::yvel].nGrowVect());

    // Creating local perturbation box array
    BoxArray m_pb_ba = turbPert.pb_ba[lev];
    Real* m_pb_um = turbPert.get_pb_um();

    // Compute u average within box union
    auto m_ixtype = xvel_data.boxArray().ixType();
    for (MFIter mfi(xvel_data, TileNoZ()) ; mfi.isValid(); ++mfi) {
        //const Box &vbx = mfi.tilebox();
        const Box &vbx = mfi.validbox();

        for (int i = 0; i < m_pb_ba.size(); i++) {
            Box pbx = convert(m_pb_ba[i], m_ixtype);
            Box ubx = pbx & vbx;
            //Box ubx = vbx & pbx; // Did not change output of ubx

            if (ubx.ok()) {
                Gpu::DeviceVector<Real> tmp_d(2,0.);
                Real* tmp = tmp_d.data();
                const Array4<const Real> & xvel_arry = xvel_data.array(mfi);

                // Operating over box union
                #ifdef USE_VOLUME_AVERAGE
                int npts = ubx.numPts();
                ParallelFor(Gpu::KernelInfo().setReduction(true), ubx, [=]
                AMREX_GPU_DEVICE(int i, int j, int k, Gpu::Handler const& handler) noexcept {
                    Gpu::deviceReduceSum(&tmp[0], xvel_arry(i,j,k), handler);
                });
                m_pb_um[i] = tmp[0] / (Real) npts;
                #endif

                #ifdef USE_SLAB_AVERAGE
                Box ubxSlab_lo = ubx.makeSlab(2,ubx.smallEnd(2));
                Box ubxSlab_hi = ubx.makeSlab(2,pbx.bigEnd(2));   //XXX pbx override
                int npts_lo = ubxSlab_lo.numPts();
                int npts_hi = ubxSlab_hi.numPts();

                Print() << "\nvbx: " << vbx
                        << "\npbx: " << pbx 
                        << "\nubx: " << ubx
                        << "\nsmallEnd: " << ubx.smallEnd() << " bigEnd: " << ubx.bigEnd()
                        << "\nubxSlab_lo[" << i << "]= " << ubxSlab_lo 
                        << "\nubxSlab_hi[" << i << "]= " << ubxSlab_hi << "\n";

                // Average in the low slab
                ParallelFor(Gpu::KernelInfo().setReduction(true), ubxSlab_lo, [=]
                AMREX_GPU_DEVICE(int i, int j, int k, Gpu::Handler const& handler) noexcept {
                    Gpu::deviceReduceSum(&tmp[0], xvel_arry(i,j,k), handler);
                });

                // Average in the high slab
                ParallelFor(Gpu::KernelInfo().setReduction(true), ubxSlab_hi, [=]
                AMREX_GPU_DEVICE(int i, int j, int k, Gpu::Handler const& handler) noexcept {
                    Gpu::deviceReduceSum(&tmp[1], xvel_arry(i,j,k), handler);
                });

                // Average the sum between top and bottom
                m_pb_um[i] = 0.5*(tmp[0] / (Real) npts_lo + tmp[1] / (Real) npts_hi); 
                #endif
            } // if
        } // for
    } // MFIter

    #ifdef DEBUG_PERTBOX_MSG
    turbPert.debug();
    #endif

    // Computing perturbation update time
    turbPert.calc_tpi_update(lev, dt);

    Print() << "Turbulent perturbation update time and amplitude initialized\n";
}

// Calculate the perturbation region amplitude.
// This function heavily emmulates the ERF::init_custom ()
void
ERF::TurbPert_amplitude (int lev, TurbulentPerturbation& turbPert)
{
    auto& lev_new = vars_new[lev];

    MultiFab r_hse(base_state[lev], make_alias, 0, 1); // r_0 is first  component
    MultiFab p_hse(base_state[lev], make_alias, 1, 1); // p_0 is second component

    MultiFab cons_pert(lev_new[Vars::cons].boxArray(), lev_new[Vars::cons].DistributionMap(),
                       lev_new[Vars::cons].nComp()   , lev_new[Vars::cons].nGrow());
    MultiFab xvel_pert(lev_new[Vars::xvel].boxArray(), lev_new[Vars::xvel].DistributionMap(), 1, lev_new[Vars::xvel].nGrowVect());
    MultiFab yvel_pert(lev_new[Vars::yvel].boxArray(), lev_new[Vars::yvel].DistributionMap(), 1, lev_new[Vars::yvel].nGrowVect());
    MultiFab zvel_pert(lev_new[Vars::zvel].boxArray(), lev_new[Vars::zvel].DistributionMap(), 1, lev_new[Vars::zvel].nGrowVect());

    // Only storing for conserved quantity for now 
    auto m_ixtype = cons_pert.boxArray().ixType();

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
        const Box &bx  = mfi.tilebox(); // Note: tilebox() = validbox() when tiling is off
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

        turbPert.apply_tpi(lev, bx, RhoTheta_comp, m_ixtype, cons_pert_arr); // Using boxArray
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
