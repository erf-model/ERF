#include <ERF_EOS.H>
#include <ERF.H>

using namespace amrex;

/**
 * Function that calls estTimeStep for each level
 *
 */
void
ERF::ComputeDt (int step)
{
    Vector<Real> dt_tmp(finest_level+1);

    for (int lev = 0; lev <= finest_level; ++lev)
    {
        dt_tmp[lev] = estTimeStep(lev, dt_mri_ratio[lev]);
    }

    ParallelDescriptor::ReduceRealMin(&dt_tmp[0], dt_tmp.size());

    Real dt_0 = dt_tmp[0];
    int n_factor = 1;
    for (int lev = 0; lev <= finest_level; ++lev) {
        dt_tmp[lev] = amrex::min(dt_tmp[lev], change_max*dt[lev]);
        n_factor *= nsubsteps[lev];
        dt_0 = amrex::min(dt_0, n_factor*dt_tmp[lev]);
        if (step == 0){
            dt_0 *= init_shrink;
            if (verbose && init_shrink != 1.0) {
                Print() << "Timestep 0: shrink initial dt at level " << lev << " by " << init_shrink << std::endl;
            }
        }
    }

    // Limit dt's by the value of stop_time.
    const Real eps = 1.e-3*dt_0;
    if (t_new[0] + dt_0 > stop_time - eps) {
        dt_0 = stop_time - t_new[0];
    }

    dt[0] = dt_0;
    for (int lev = 1; lev <= finest_level; ++lev) {
        dt[lev] = dt[lev-1] / nsubsteps[lev];
    }
}

/**
 * Function that calls estTimeStep for each level
 *
 * @param[in] level level of refinement (coarsest level i 0)
 * @param[out] dt_fast_ratio ratio of slow to fast time step
 */
Real
ERF::estTimeStep (int level, long& dt_fast_ratio) const
{
    BL_PROFILE("ERF::estTimeStep()");

    Real estdt_comp = 1.e20;
    Real estdt_lowM = 1.e20;

    auto const dxinv = geom[level].InvCellSizeArray();
    auto const dzinv = 1.0 / dz_min;

    MultiFab const& S_new = vars_new[level][Vars::cons];

    MultiFab ccvel(grids[level],dmap[level],3,0);

    average_face_to_cellcenter(ccvel,0,
                               Array<const MultiFab*,3>{&vars_new[level][Vars::xvel],
                                                        &vars_new[level][Vars::yvel],
                                                        &vars_new[level][Vars::zvel]});

    int l_no_substepping = solverChoice.no_substepping;
    int l_incompressible = solverChoice.incompressible[level];

#ifdef ERF_USE_EB
    EBFArrayBoxFactory ebfact = EBFactory(level);
    const MultiFab& detJ = ebfact.getVolFrac();
#endif

#ifdef ERF_USE_EB
    Real estdt_comp_inv = ReduceMax(S_new, ccvel, detJ, 0,
       [=] AMREX_GPU_HOST_DEVICE (Box const& b,
                                  Array4<Real const> const& s,
                                  Array4<Real const> const& vf,
                                  Array4<Real const> const& u) -> Real
#else
    Real estdt_comp_inv = ReduceMax(S_new, ccvel, 0,
       [=] AMREX_GPU_HOST_DEVICE (Box const& b,
                                  Array4<Real const> const& s,
                                  Array4<Real const> const& u) -> Real
#endif
       {
           Real new_comp_dt = -1.e100;
           amrex::Loop(b, [=,&new_comp_dt] (int i, int j, int k) noexcept
           {
#ifdef ERF_USE_EB
               if (vf(i,j,k) > 0.)
#endif
               {
                   const Real rho      = s(i, j, k, Rho_comp);
                   const Real rhotheta = s(i, j, k, RhoTheta_comp);

                   // NOTE: even when moisture is present,
                   //       we only use the partial pressure of the dry air
                   //       to compute the soundspeed
                   Real pressure = getPgivenRTh(rhotheta);
                   Real c = std::sqrt(Gamma * pressure / rho);

                   // If we are not doing the acoustic substepping, then the z-direction contributes
                   //    to the computation of the time step
                   if (l_no_substepping) {
                       new_comp_dt = amrex::max(((amrex::Math::abs(u(i,j,k,0))+c)*dxinv[0]),
                                                ((amrex::Math::abs(u(i,j,k,1))+c)*dxinv[1]),
                                                ((amrex::Math::abs(u(i,j,k,2))+c)*dzinv   ), new_comp_dt);

                   // If we are     doing the acoustic substepping, then the z-direction does not contribute
                   //    to the computation of the time step
                   } else {
                       new_comp_dt = amrex::max(((amrex::Math::abs(u(i,j,k,0))+c)*dxinv[0]),
                                                ((amrex::Math::abs(u(i,j,k,1))+c)*dxinv[1]), new_comp_dt);
                   }
               }
           });
           return new_comp_dt;
       });

    ParallelDescriptor::ReduceRealMax(estdt_comp_inv);
    estdt_comp = cfl / estdt_comp_inv;

     Real estdt_lowM_inv = ReduceMax(ccvel, 0,
       [=] AMREX_GPU_HOST_DEVICE (Box const& b,
                                  Array4<Real const> const& u) -> Real
       {
           Real new_lm_dt = -1.e100;
           Loop(b, [=,&new_lm_dt] (int i, int j, int k) noexcept
           {
               new_lm_dt = amrex::max(((amrex::Math::abs(u(i,j,k,0)))*dxinv[0]),
                                      ((amrex::Math::abs(u(i,j,k,1)))*dxinv[1]),
                                      ((amrex::Math::abs(u(i,j,k,2)))*dxinv[2]), new_lm_dt);
           });
           return new_lm_dt;
       });

     ParallelDescriptor::ReduceRealMax(estdt_lowM_inv);
     if (estdt_lowM_inv > 0.0_rt)
         estdt_lowM = cfl / estdt_lowM_inv;

     if (verbose) {
         if (fixed_dt[level] <= 0.0) {
             Print() << "Using cfl = " << cfl << std::endl;
             Print() << "Compressible dt at level " << level << ":  " << estdt_comp << std::endl;
             if (estdt_lowM_inv > 0.0_rt) {
                 Print() << "Incompressible dt at level " << level << ":  " << estdt_lowM << std::endl;
             } else {
                 Print() << "Incompressible dt at level " << level << ": undefined " << std::endl;
             }
         }

         if (fixed_dt[level] > 0.0) {
             Print() << "Based on cfl of 1.0 " << std::endl;
             Print() << "Compressible dt at level " << level << " would be:  " << estdt_comp/cfl << std::endl;
             if (estdt_lowM_inv > 0.0_rt) {
                 Print() << "Incompressible dt at level " << level << " would be:  " << estdt_lowM/cfl << std::endl;
             } else {
                 Print() << "Incompressible dt at level " << level << " would be undefined " << std::endl;
             }
             Print() << "Fixed dt at level " << level << "       is:  " << fixed_dt[level] << std::endl;
             if (fixed_fast_dt[level] > 0.0) {
                 Print() << "Fixed fast dt at level " << level << "       is:  " << fixed_fast_dt[level] << std::endl;
             }
         }
     }

     if (!l_no_substepping) {
         if (fixed_dt[level] > 0. && fixed_fast_dt[level] > 0.) {
             dt_fast_ratio = static_cast<long>( fixed_dt[level] / fixed_fast_dt[level] );
         } else if (fixed_dt[level] > 0.) {
             dt_fast_ratio = static_cast<long>( std::ceil((fixed_dt[level]/estdt_comp)) );
         } else {
             dt_fast_ratio = (estdt_lowM_inv > 0.0) ? static_cast<long>( std::ceil((estdt_lowM/estdt_comp)) ) : 1;
         }

         // Force time step ratio to be an even value
         if (solverChoice.force_stage1_single_substep) {
             if ( dt_fast_ratio%2 != 0) dt_fast_ratio += 1;
         } else {
             if ( dt_fast_ratio%6 != 0) {
                 Print() << "mri_dt_ratio = " << dt_fast_ratio
                         << " not divisible by 6 for N/3 substeps in stage 1" << std::endl;
                 dt_fast_ratio = static_cast<int>(std::ceil(dt_fast_ratio/6.0) * 6);
             }
         }

         if (verbose) {
             Print() << "smallest even ratio is: " << dt_fast_ratio << std::endl;
         }
     } // if substepping

     if (fixed_dt[level] > 0.0) {
         return fixed_dt[level];
     } else {
         // Incompressible (substepping is not allowed)
         if (l_incompressible) {
             return estdt_lowM;

         // Compressible with or without substepping
         } else {
             return estdt_comp;
         }
     }
}
