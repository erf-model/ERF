#ifndef ERF_SRI_H
#define ERF_SRI_H

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_ParmParse.H>
#include <AMReX_IntegratorBase.H>
#include <TimeIntegration.H>
#include <functional>

namespace SRI {
    enum struct ButcherTableauTypes {
        User = 0,
        ForwardEuler,
        Trapezoid,
        SSPRK3,
        RK4,
        NumTypes
    };
}

template<class T>
class SRIIntegrator
{
private:
   /**
    * \brief rhs is the right-hand-side function the integrator will use.
    */
    std::function<void(T&, const T&, T&, const amrex::Real, const amrex::Real, int)> rhs;

   /**
    * \brief Integrator timestep size (Real)
    */
    amrex::Real timestep;

   /**
    * \brief Should we allocate space for and update fluxes
    */
    int use_fluxes;

   /**
    * \brief The post_update function is called by the integrator on state data before using it to evaluate a right-hand side.
    */
    std::function<void (T&, amrex::Real, amrex::Real, int, int)> post_update;

    SRI::ButcherTableauTypes butcher_type;
    int number_nodes;
    amrex::Vector<std::unique_ptr<T> > F_nodes;
    amrex::Vector<amrex::Vector<amrex::Real> > tableau;
    amrex::Vector<amrex::Gpu::DeviceVector<amrex::Real> > tableau_d;
    amrex::Vector<amrex::Real> weights;
    amrex::Gpu::DeviceVector<amrex::Real> weights_d;
    amrex::Vector<amrex::Real> nodes;

    void initialize_preset_tableau (SRI::ButcherTableauTypes tableau_type)
    {
        butcher_type = tableau_type;
        switch (tableau_type)
        {
            case SRI::ButcherTableauTypes::ForwardEuler:
                nodes = {0.0};
                tableau = {{0.0}};
                weights = {1.0};
                break;
            case SRI::ButcherTableauTypes::Trapezoid:
                nodes = {0.0,
                        1.0};
                tableau = {{0.0},
                        {1.0, 0.0}};
                weights = {0.5, 0.5};
                break;
            case SRI::ButcherTableauTypes::SSPRK3:
                nodes = {0.0,
                        1.0,
                        0.5};
                tableau = {{0.0},
                        {1.0, 0.0},
                        {0.25, 0.25, 0.0}};
                weights = {1./6., 1./6., 2./3.};
                break;
            case SRI::ButcherTableauTypes::RK4:
                nodes = {0.0,
                        0.5,
                        0.5,
                        1.0};
                tableau = {{0.0},
                        {0.5, 0.0},
                        {0.0, 0.5, 0.0},
                        {0.0, 0.0, 1.0, 0.0}};
                weights = {1./6., 1./3., 1./3., 1./6.};
                break;
            default:
                amrex::Error("Invalid RK Integrator tableau type integration.rk.type");
                break;
        }

        number_nodes = weights.size();
        weights_d.resize(weights.size());
        amrex::Gpu::copy(amrex::Gpu::hostToDevice, weights.begin(), weights.end(),
                         weights_d.begin());

        tableau_d.resize(tableau.size());
        for (int i = 0; i < number_nodes; ++i) {
            tableau_d[i].resize(tableau[i].size());
            amrex::Gpu::copy(amrex::Gpu::hostToDevice, tableau[i].begin(), tableau[i].end(),
                            tableau_d[i].begin());
        }
    }

    void initialize_data (const T& S_data)
    {
        // Create data for stage RHS
        const bool include_ghost = true;
        for (int i = 0; i < number_nodes; ++i)
        {
            amrex::IntegratorOps<T>::CreateLike(F_nodes, S_data, include_ghost);
        }
    }

public:
    SRIIntegrator () {}

    SRIIntegrator (const T& S_data)
    {
        initialize(S_data);
    }

    void initialize (const T& S_data)
    {
        amrex::ParmParse pp("integration.rk");

        // Read an integrator type, if not recognized, then read weights/nodes/butcher tableau
        // By default, use SSPRK3
        int _tableau_type = static_cast<int>(SRI::ButcherTableauTypes::SSPRK3);
        pp.query("type", _tableau_type);

        initialize_preset_tableau(static_cast<SRI::ButcherTableauTypes>(_tableau_type));
        initialize_data(S_data);
    }

    virtual ~SRIIntegrator () {}

    void setUseFluxes(int _use_fluxes)
    {
        use_fluxes = _use_fluxes;
    }

    void set_rhs (std::function<void(T&, const T&, T&, const amrex::Real, const amrex::Real, int)> F)
    {
        rhs = F;
    }

    void set_post_update (std::function<void (T&, amrex::Real, amrex::Real, int, int)> F)
    {
        post_update = F;
    }

    std::function<void (T&, amrex::Real, amrex::Real, int, int)> get_post_update ()
    {
        return post_update;
    }

    std::function<void(T&, const T&, T&, const amrex::Real, int)> get_rhs ()
    {
        return rhs;
    }

    amrex::Real advance (T& S_old, T& S_new, amrex::Real time, const amrex::Real time_step)
    {
        using namespace amrex;

        // Assume before advance() that S_old is valid data at the current time ("time" argument)
        // And that both S_old and S_new contain ghost cells for evaluating a stencil based RHS

        // **********************************************
        // RK Integration
        // **********************************************

        // We copy all variables when we call Copy in the following -- the velocities
        //    will be over-written in the rhs on all valid faces but we
        //    use this copy to fill in the ghost locations which will
        //    be needed for metric terms

        timestep = time_step;

        int nav = Cons::NumVars;
        const amrex::GpuArray<int, 7> scomp_all = {0,0,0,0,0,0,0};
        const amrex::GpuArray<int, 7> ncomp_all = {nav,1,1,1,nav,nav,nav};

        amrex::Vector<amrex::MultiFab> S_scratch;
        Real* p_weights_d = weights_d.dataPtr();

        // Fill the RHS F_nodes at each stage
        for (int in = 0; in < number_nodes; ++in)
        {
            // Get current stage time, t = t_old + h * Ci
            amrex::Real dt_stage   = time_step * nodes[in];
            amrex::Real stage_time = time + dt_stage;

            // Fill S_new with the solution value for evaluating F at the current stage
            // Copy S_new = S_old then Saxpy across the butcher table row:
            // S_new += h * Aij * Fj

            Real* tableau_in_p = tableau_d[in].dataPtr();

    #ifdef _OPENMP
    #pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
    #endif
            {
                for ( MFIter mfi(S_new[IntVar::cons],TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                    const Box& bx = mfi.tilebox();
                    const Box& tbx = mfi.nodaltilebox(0);
                    const Box& tby = mfi.nodaltilebox(1);
                    const Box& tbz = mfi.nodaltilebox(2);

                    const Array4<Real>& snew_cons = S_new[IntVar::cons].array(mfi);
                    const Array4<Real>& snew_xmom = S_new[IntVar::xmom].array(mfi);
                    const Array4<Real>& snew_ymom = S_new[IntVar::ymom].array(mfi);
                    const Array4<Real>& snew_zmom = S_new[IntVar::zmom].array(mfi);

                    const Array4<Real>& sold_cons = S_old[IntVar::cons].array(mfi);
                    const Array4<Real>& sold_xmom = S_old[IntVar::xmom].array(mfi);
                    const Array4<Real>& sold_ymom = S_old[IntVar::ymom].array(mfi);
                    const Array4<Real>& sold_zmom = S_old[IntVar::zmom].array(mfi);

                    Vector<Array4<Real > > frhs_arr;
                    Gpu::DeviceVector<Array4<Real > > frhs_arr_d;
                    for (auto& f : F_nodes) {
                        for (auto& ff : *f) {
                            frhs_arr.push_back(ff.array(mfi));
                        }
                    }
                    frhs_arr_d.resize(frhs_arr.size());
                    amrex::Gpu::copy(amrex::Gpu::hostToDevice, frhs_arr.begin(), frhs_arr.end(),
                                    frhs_arr_d.begin());
                    Array4<Real>* p_frhs_arr_d = frhs_arr_d.dataPtr();
                    const int fstride = S_new.size();
                    const int icons = 0;
                    const int ixmom = 1;
                    const int iymom = 2;
                    const int izmom = 3;

                    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                        for (int n = scomp_all[icons]; n < scomp_all[icons] + ncomp_all[icons]; ++n) {
                            snew_cons(i,j,k,n) = sold_cons(i,j,k,n);
                            for (int fn = 0; fn < in; ++fn) {
                                snew_cons(i,j,k,n) += time_step * tableau_in_p[fn] * p_frhs_arr_d[fn*fstride+icons](i,j,k,n);
                            }
                        }
                    });

                    ParallelFor(tbx, tby, tbz,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                        for (int n = scomp_all[ixmom]; n < scomp_all[ixmom] + ncomp_all[ixmom]; ++n) {
                            snew_xmom(i,j,k,n) = sold_xmom(i,j,k,n);
                            for (int fn = 0; fn < in; ++fn) {
                                snew_xmom(i,j,k,n) += time_step * tableau_in_p[fn] * p_frhs_arr_d[fn*fstride+ixmom](i,j,k,n);
                            }
                        }
                    },
                    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                        for (int n = scomp_all[iymom]; n < scomp_all[iymom] + ncomp_all[iymom]; ++n) {
                            snew_ymom(i,j,k,n) = sold_ymom(i,j,k,n);
                            for (int fn = 0; fn < in; ++fn) {
                                snew_ymom(i,j,k,n) += time_step * tableau_in_p[fn] * p_frhs_arr_d[fn*fstride+iymom](i,j,k,n);
                            }
                        }
                    },
                    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                        for (int n = scomp_all[izmom]; n < scomp_all[izmom] + ncomp_all[izmom]; ++n) {
                            snew_zmom(i,j,k,n) = sold_zmom(i,j,k,n);
                            for (int fn = 0; fn < in; ++fn) {
                                snew_zmom(i,j,k,n) += time_step * tableau_in_p[fn] * p_frhs_arr_d[fn*fstride+izmom](i,j,k,n);
                            }
                        }
                    });

                    if (use_fluxes) {
                        const int ixflx = 4;
                        const int iyflx = 5;
                        const int izflx = 6;

                        const Array4<Real>& snew_xflx = S_new[IntVar::xflux].array(mfi);
                        const Array4<Real>& snew_yflx = S_new[IntVar::yflux].array(mfi);
                        const Array4<Real>& snew_zflx = S_new[IntVar::zflux].array(mfi);
                        const Array4<Real>& sold_xflx = S_old[IntVar::xflux].array(mfi);
                        const Array4<Real>& sold_yflx = S_old[IntVar::yflux].array(mfi);
                        const Array4<Real>& sold_zflx = S_old[IntVar::zflux].array(mfi);

                        ParallelFor(tbx, tby, tbz,
                        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                            for (int n = scomp_all[ixflx]; n < scomp_all[ixflx] + ncomp_all[ixflx]; ++n) {
                                snew_xflx(i,j,k,n) = sold_xflx(i,j,k,n);
                                for (int fn = 0; fn < in; ++fn) {
                                    snew_xflx(i,j,k,n) += time_step * tableau_in_p[fn] * p_frhs_arr_d[fn*fstride+ixflx](i,j,k,n);
                                }
                            }
                        },
                            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                            for (int n = scomp_all[iyflx]; n < scomp_all[iyflx] + ncomp_all[iyflx]; ++n) {
                                snew_yflx(i,j,k,n) = sold_yflx(i,j,k,n);
                                for (int fn = 0; fn < in; ++fn) {
                                    snew_yflx(i,j,k,n) += time_step * tableau_in_p[fn] * p_frhs_arr_d[fn*fstride+iyflx](i,j,k,n);
                                }
                            }
                        },
                        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                            for (int n = scomp_all[izflx]; n < scomp_all[izflx] + ncomp_all[izflx]; ++n) {
                                snew_zflx(i,j,k,n) = sold_zflx(i,j,k,n);
                                for (int fn = 0; fn < in; ++fn) {
                                snew_zflx(i,j,k,n) += time_step * tableau_in_p[fn] * p_frhs_arr_d[fn*fstride+izflx](i,j,k,n);
                                }
                            }
                        });
                    }

                }
            }

            // Call the post-update hook for the stage state value
            post_update(S_new, stage_time, dt_stage, S_new[IntVar::cons].nGrow(), S_new[IntVar::xmom].nGrow());

            // Fill F[i], the RHS at the current stage
            // F[i] = RHS(y, t) at y = stage_value, t = stage_time
            rhs(*F_nodes[in], S_new, S_scratch, time, stage_time, RHSVar::all);
        }

        // Fill new State, starting with S_new = S_old.
        // Then Saxpy S_new += h * Wi * Fi for integration weights Wi
        // The following MFIter simply fuses these operations into a single MFIter loop.
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        {
            for ( MFIter mfi(S_new[IntVar::cons],TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                const Box& bx = mfi.tilebox();
                const Box& tbx = mfi.nodaltilebox(0);
                const Box& tby = mfi.nodaltilebox(1);
                const Box& tbz = mfi.nodaltilebox(2);

                const Array4<Real>& snew_cons = S_new[IntVar::cons].array(mfi);
                const Array4<Real>& snew_xmom = S_new[IntVar::xmom].array(mfi);
                const Array4<Real>& snew_ymom = S_new[IntVar::ymom].array(mfi);
                const Array4<Real>& snew_zmom = S_new[IntVar::zmom].array(mfi);

                const Array4<Real>& sold_cons = S_old[IntVar::cons].array(mfi);
                const Array4<Real>& sold_xmom = S_old[IntVar::xmom].array(mfi);
                const Array4<Real>& sold_ymom = S_old[IntVar::ymom].array(mfi);
                const Array4<Real>& sold_zmom = S_old[IntVar::zmom].array(mfi);

                Vector<Array4<Real > > frhs_arr;
                Gpu::DeviceVector<Array4<Real > > frhs_arr_d;
                for (auto& f : F_nodes) {
                    for (auto& ff : *f) {
                        frhs_arr.push_back(ff.array(mfi));
                    }
                }
                frhs_arr_d.resize(frhs_arr.size());
                amrex::Gpu::copy(amrex::Gpu::hostToDevice, frhs_arr.begin(), frhs_arr.end(),
                                frhs_arr_d.begin());
                Array4<Real>* p_frhs_arr_d = frhs_arr_d.dataPtr();
                const int fstride = S_new.size();
                const int icons = 0;
                const int ixmom = 1;
                const int iymom = 2;
                const int izmom = 3;
                const int nnodes = number_nodes;

                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                    for (int n = scomp_all[icons]; n < scomp_all[icons] + ncomp_all[icons]; ++n) {
                        snew_cons(i,j,k,n) = sold_cons(i,j,k,n);
                        for (int fn = 0; fn < nnodes; ++fn) {
                            snew_cons(i,j,k,n) += time_step * p_weights_d[fn] * p_frhs_arr_d[fn*fstride+icons](i,j,k,n);
                        }
                    }
                });

                ParallelFor(tbx, tby, tbz,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                    for (int n = scomp_all[ixmom]; n < scomp_all[ixmom] + ncomp_all[ixmom]; ++n) {
                        snew_xmom(i,j,k,n) = sold_xmom(i,j,k,n);
                        for (int fn = 0; fn < nnodes; ++fn) {
                            snew_xmom(i,j,k,n) += time_step * p_weights_d[fn] * p_frhs_arr_d[fn*fstride+ixmom](i,j,k,n);
                        }
                    }
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                    for (int n = scomp_all[iymom]; n < scomp_all[iymom] + ncomp_all[iymom]; ++n) {
                        snew_ymom(i,j,k,n) = sold_ymom(i,j,k,n);
                        for (int fn = 0; fn < nnodes; ++fn) {
                            snew_ymom(i,j,k,n) += time_step * p_weights_d[fn] * p_frhs_arr_d[fn*fstride+iymom](i,j,k,n);
                        }
                    }
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                    for (int n = scomp_all[izmom]; n < scomp_all[izmom] + ncomp_all[izmom]; ++n) {
                        snew_zmom(i,j,k,n) = sold_zmom(i,j,k,n);
                        for (int fn = 0; fn < nnodes; ++fn) {
                            snew_zmom(i,j,k,n) += time_step * p_weights_d[fn] * p_frhs_arr_d[fn*fstride+izmom](i,j,k,n);
                        }
                    }
                });

                if (use_fluxes) {
                    const int ixflx = 4;
                    const int iyflx = 5;
                    const int izflx = 6;

                    const Array4<Real>& snew_xflx = S_new[IntVar::xflux].array(mfi);
                    const Array4<Real>& snew_yflx = S_new[IntVar::yflux].array(mfi);
                    const Array4<Real>& snew_zflx = S_new[IntVar::zflux].array(mfi);
                    const Array4<Real>& sold_xflx = S_old[IntVar::xflux].array(mfi);
                    const Array4<Real>& sold_yflx = S_old[IntVar::yflux].array(mfi);
                    const Array4<Real>& sold_zflx = S_old[IntVar::zflux].array(mfi);

                    ParallelFor(tbx, tby, tbz,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                        for (int n = scomp_all[ixflx]; n < scomp_all[ixflx] + ncomp_all[ixflx]; ++n) {
                            snew_xflx(i,j,k,n) = sold_xflx(i,j,k,n);
                            for (int fn = 0; fn < nnodes; ++fn) {
                                snew_xflx(i,j,k,n) += time_step * p_weights_d[fn] * p_frhs_arr_d[fn*fstride+ixflx](i,j,k,n);
                            }
                        }
                    },
                    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                            for (int n = scomp_all[iyflx]; n < scomp_all[iyflx] + ncomp_all[iyflx]; ++n) {
                            snew_yflx(i,j,k,n) = sold_yflx(i,j,k,n);
                            for (int fn = 0; fn < nnodes; ++fn) {
                                snew_yflx(i,j,k,n) += time_step * p_weights_d[fn] * p_frhs_arr_d[fn*fstride+iyflx](i,j,k,n);
                            }
                        }
                    },
                    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                        for (int n = scomp_all[izflx]; n < scomp_all[izflx] + ncomp_all[izflx]; ++n) {
                            snew_zflx(i,j,k,n) = sold_zflx(i,j,k,n);
                            for (int fn = 0; fn < nnodes; ++fn) {
                                snew_zflx(i,j,k,n) += time_step * p_weights_d[fn] * p_frhs_arr_d[fn*fstride+izflx](i,j,k,n);
                            }
                        }
                    });
                }
            }
        }

        // Call the post-update hook for S_new
        post_update(S_new, time + time_step, time_step, S_new[IntVar::cons].nGrow(), S_new[IntVar::xmom].nGrow());

        // Return timestep
        return time_step;
    }

    void time_interpolate (const T& /* S_new */, const T& /*S_old*/, amrex::Real /*timestep_fraction*/, T& /*data*/) {}

    void map_data (std::function<void(T&)> Map)
    {
        for (auto& F : F_nodes) {
            Map(*F);
        }
    }

};
#endif
