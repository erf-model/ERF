        // *****************************************************************************
        // Limit the slow scalars once the whole stage update is complete: after the
        //     explicit update in erf_slow_rhs_post, the implicit vertical diffusion
        //     in ERF_ImplicitPost.H and (with EB) state redistribution.  Doing this in
        //     one place keeps explicit and implicit diffusion consistent.
        //  - rho*KE is floored at rho*tke_floor (or machine epsilon if tke_floor <= 0)
        //  - rho*Q  is floored at zero
        // *****************************************************************************
        {
            BL_PROFILE("slow_rhs_post_limit");

            const TurbChoice& l_tc = solverChoice.turbChoice[level];

            const bool l_limit_KE = l_tc.use_tke;
            const bool l_limit_Q  = (solverChoice.moisture_type != MoistureType::None);
            const bool l_sbm      = (solverChoice.moisture_type == MoistureType::SBM);

            // The Dirichlet wall value of k (first cell above the wall) is pinned by
            // erf_slow_rhs_post and the implicit solve; leave it untouched here.
            const bool l_dirichlet_k = ( l_tc.rans_type == RANSType::kEqn && l_tc.dirichlet_k &&
                                         (m_SurfaceLayer[Orientation(Direction::z, Orientation::low)] != nullptr) );
            const int  l_klo         = fine_geom.Domain().smallEnd(2);

            const Real l_tke_floor = l_tc.tke_floor;
            const Real l_eps       = std::numeric_limits<Real>::epsilon();

            // Every state component from RhoQ1_comp to the end of the state.  With SBM
            // only the host-owned qv is limited; the projected liquid lanes are owned by
            // the spectral solver.
            const int l_nq = (l_sbm) ? 1 : S_new[IntVars::cons].nComp() - RhoQ1_comp;

            if (l_limit_Q) {
                erf_sbm::require_host_write_range_allowed(
                    l_sbm, RhoQ1_comp, l_nq, erf_sbm::HostWritePath::Positivity);
            }

            if (l_limit_KE || l_limit_Q) {
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
                for ( MFIter mfi(S_new[IntVars::cons],TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    const Box& tbx = mfi.tilebox();
                    const Array4<Real>& cons = S_new[IntVars::cons].array(mfi);

                    if (l_limit_KE) {
                        const bool skip_klo = l_dirichlet_k && (tbx.smallEnd(2) <= l_klo);
                        ParallelFor(tbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                        {
                            if (skip_klo && k == l_klo) { return; }
                            const Real ke_floor = (l_tke_floor > zero) ? cons(i,j,k,Rho_comp) * l_tke_floor : l_eps;
                            cons(i,j,k,RhoKE_comp) = amrex::max(cons(i,j,k,RhoKE_comp), ke_floor);
                        });
                    }

                    if (l_limit_Q) {
                        ParallelFor(tbx, l_nq, [=,zero_d=zero] AMREX_GPU_DEVICE (int i, int j, int k, int nn) noexcept
                        {
                            const int n = RhoQ1_comp + nn;
                            cons(i,j,k,n) = amrex::max(cons(i,j,k,n), zero_d);
                        });
                    }
                } // mfi
            } // limit KE or Q
        } // profile
