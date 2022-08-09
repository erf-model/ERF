#ifndef ERF_MRI_H
#define ERF_MRI_H
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_ParmParse.H>
#include <AMReX_IntegratorBase.H>
#include <TimeIntegration.H>
#include <functional>

template<class T>
class MRISplitIntegrator : public amrex::IntegratorBase<T>
{
private:
   /**
    * \brief rhs is the right-hand-side function the integrator will use.
    */
    std::function<void(T&,         const T&,     const amrex::Real, int)> rhs;
    std::function<void(T&,         const T&, T&, const amrex::Real, int)> slow_rhs;
    std::function<void(T&, T&, T&, const T&, T&, const amrex::Real, const amrex::Real)> fast_rhs;

   /**
    * \brief Integrator timestep size (Real)
    */
    amrex::Real timestep;

   /**
    * \brief The ratio of slow timestep size / fast timestep size (int)
    */
    int slow_fast_timestep_ratio = 0;

   /**
    * \brief The post_update function is called by the integrator on state data before using it to evaluate a right-hand side.
    */
    std::function<void (T&, amrex::Real, int, int)> post_update;

    std::function<void (T&, amrex::Real, int, int)> post_substep;


    amrex::Vector<std::unique_ptr<T> > T_store;
    T* S_sum;
    T* S_scratch;
    T* F_slow;
    T* F_pert;

    void initialize_data (const T& S_data)
    {
        // TODO: We can optimize memory by making the cell-centered part of S_sum, S_scratch and F_pert
        //       have only 2 components, not Cons::NumVars components
        const bool include_ghost = true;
        amrex::IntegratorOps<T>::CreateLike(T_store, S_data, include_ghost);
        S_sum = T_store[0].get();
        amrex::IntegratorOps<T>::CreateLike(T_store, S_data, include_ghost);
        S_scratch = T_store[1].get();
        amrex::IntegratorOps<T>::CreateLike(T_store, S_data, include_ghost);
        F_slow = T_store[2].get();
        amrex::IntegratorOps<T>::CreateLike(T_store, S_data, include_ghost);
        F_pert = T_store[3].get();
    }

public:
    MRISplitIntegrator () {}

    MRISplitIntegrator (const T& S_data)
    {
        initialize_data(S_data);
    }

    void initialize (const T& S_data)
    {
        initialize_data(S_data);
    }

    virtual ~MRISplitIntegrator () {}

    void set_rhs (std::function<void(T&, const T&, const amrex::Real, int)> F)
    {
        rhs = F;
    }

    void set_slow_rhs (std::function<void(T&, const T&, T&, const amrex::Real, int)> F)
    {
        slow_rhs = F;
    }

    void set_fast_rhs (std::function<void(T&, T&, T&, const T&, T&, const amrex::Real, const amrex::Real)> F)
    {
        fast_rhs = F;
    }

    std::function<void(T&, T&, T&, const T&, const amrex::Real, const amrex::Real)> get_fast_rhs ()
    {
        return fast_rhs;
    }

    void set_slow_fast_timestep_ratio (const int timestep_ratio = 1)
    {
        slow_fast_timestep_ratio = timestep_ratio;
    }

    int get_slow_fast_timestep_ratio ()
    {
        return slow_fast_timestep_ratio;
    }

    void set_post_update (std::function<void (T&, amrex::Real, int, int)> F)
    {
        post_update = F;
    }

    void set_post_substep (std::function<void (T&, amrex::Real, int, int)> F)
    {
        post_substep = F;
    }

    std::function<void (T&, amrex::Real, int, int)> get_post_update ()
    {
        return post_update;
    }

    std::function<void(T&, const T&, const amrex::Real, int)> get_rhs ()
    {
        return rhs;
    }

    std::function<void (T&, amrex::Real, int, int)> get_post_substep ()
    {
        return post_substep;
    }

    amrex::Real advance (T& S_old, T& S_new, amrex::Real time, const amrex::Real time_step)
    {
        using namespace amrex;

        timestep = time_step;

        const int substep_ratio = get_slow_fast_timestep_ratio();

        AMREX_ALWAYS_ASSERT(substep_ratio > 1 && substep_ratio % 2 == 0);

        const amrex::Real sub_timestep = timestep / substep_ratio;

        // Assume before advance() that S_old is valid data at the current time ("time" argument)
        // And that if data is a MultiFab, both S_old and S_new contain ghost cells for evaluating a stencil based RHS
        // We need this from S_old. This is convenient for S_new to have so we can use it
        // as scratch space for stage values without creating a new scratch MultiFab with ghost cells.

        // NOTE: In the following, we use S_new to hold S*, S**, and finally, S^(n+1) at the new time
        // DEFINITIONS:
        // S_old  = S^n
        // S_sum  = S(t)
        // F_slow = F(S_stage)
        // F_pert = G(S(t)-S_stage, S_stage)

        /**********************************************/
        /* RK3 Integration with Acoustic Sub-stepping */
        /**********************************************/

        // Start with S_new (aka S_stage) holding S_old
    #ifdef _OPENMP
    #pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
    #endif
        {
            for ( MFIter mfi(S_old[IntVar::cons],TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                const Box gbx = mfi.tilebox().grow(S_old[IntVar::cons].nGrowVect());
                const Box gtbx = mfi.nodaltilebox(0).grow(S_old[IntVar::xmom].nGrowVect());
                const Box gtby = mfi.nodaltilebox(1).grow(S_old[IntVar::ymom].nGrowVect());
                const Box gtbz = mfi.nodaltilebox(2).grow(S_old[IntVar::zmom].nGrowVect());
                const Box gfbx = mfi.nodaltilebox(0).grow(S_old[IntVar::xflux].nGrowVect());
                const Box gfby = mfi.nodaltilebox(1).grow(S_old[IntVar::yflux].nGrowVect());
                const Box gfbz = mfi.nodaltilebox(2).grow(S_old[IntVar::zflux].nGrowVect());

                Vector<Array4<Real> > sold_h(IntVar::NumVars);
                Vector<Array4<Real> > snew_h(IntVar::NumVars);

                for (int i = 0; i < IntVar::NumVars; ++i) {
                    sold_h[i] = S_old[i].array(mfi);
                    snew_h[i] = S_new[i].array(mfi);
                }

                Gpu::AsyncVector<Array4<Real> > sold_d(IntVar::NumVars);
                Gpu::AsyncVector<Array4<Real> > snew_d(IntVar::NumVars);

                Gpu::copy(Gpu::hostToDevice, sold_h.begin(), sold_h.end(), sold_d.begin());
                Gpu::copy(Gpu::hostToDevice, snew_h.begin(), snew_h.end(), snew_d.begin());

                Array4<Real>* sold = sold_d.dataPtr();
                Array4<Real>* snew = snew_d.dataPtr();

                ParallelFor(gbx, Cons::NumVars,
                [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
                    snew[IntVar::cons](i,j,k,n) = sold[IntVar::cons](i,j,k,n);
                });

                ParallelFor(gtbx, gtby, gtbz,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                    snew[IntVar::xmom](i,j,k) = sold[IntVar::xmom](i,j,k);
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                    snew[IntVar::ymom](i,j,k) = sold[IntVar::ymom](i,j,k);
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                    snew[IntVar::zmom](i,j,k) = sold[IntVar::zmom](i,j,k);
                });

                ParallelFor(
                gfbx, Cons::NumVars,
                [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept {
                    snew[IntVar::xflux](i,j,k,n) = sold[IntVar::xflux](i,j,k,n);
                },
                gfby, Cons::NumVars,
                [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept {
                    snew[IntVar::yflux](i,j,k,n) = sold[IntVar::yflux](i,j,k,n);
                },
                gfbz, Cons::NumVars,
                [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept {
                    snew[IntVar::zflux](i,j,k,n) = sold[IntVar::zflux](i,j,k,n);
                });
            }
        }
        post_update(S_new, time, S_new[IntVar::cons].nGrow(), S_new[IntVar::xmom].nGrow());

        // Timestep taken by the fast integrator
        amrex::Real dtau;

        // How many timesteps taken by the fast integrator
        int nsubsteps;

        int nav = Cons::NumVars;
        const amrex::GpuArray<int, IntVar::NumVars> scomp_all = {0,0,0,0,0,0,0};
        const amrex::GpuArray<int, IntVar::NumVars> ncomp_all = {nav,1,1,1,nav,nav,nav};

        const amrex::GpuArray<int, IntVar::NumVars> scomp_fast = {0,0,0,0,0,0,0};
        const amrex::GpuArray<int, IntVar::NumVars> ncomp_fast = {2,1,1,1,2,2,2};

        int nsv = Cons::NumVars-2;
        const amrex::GpuArray<int, IntVar::NumVars> scomp_slow = {  2,0,0,0,  2,  2,  2};
        const amrex::GpuArray<int, IntVar::NumVars> ncomp_slow = {nsv,0,0,0,nsv,nsv,nsv};

        const amrex::GpuArray<int, IntVar::NumVars> scomp_rth  = {1,0,0,0,0,0,0};
        const amrex::GpuArray<int, IntVar::NumVars> ncomp_rth  = {1,0,0,0,0,0,0};

        // We copy (rho theta) and the velocities here -- the velocities
        //    will be over-written in slow_rhs on all valid faces but we
        //    use this copy to fill in the ghost locations which will
        //    be needed for metric terms
    #ifdef _OPENMP
    #pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
    #endif
        {
            for ( MFIter mfi(S_old[IntVar::cons],TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                const Box gtbx = mfi.nodaltilebox(0).grow(S_old[IntVar::xmom].nGrowVect());
                const Box gtby = mfi.nodaltilebox(1).grow(S_old[IntVar::ymom].nGrowVect());
                const Box gtbz = mfi.nodaltilebox(2).grow(S_old[IntVar::zmom].nGrowVect());

                const Array4<Real>& scrh_xmom = (*S_scratch)[IntVar::xmom].array(mfi);
                const Array4<Real>& scrh_ymom = (*S_scratch)[IntVar::ymom].array(mfi);
                const Array4<Real>& scrh_zmom = (*S_scratch)[IntVar::zmom].array(mfi);

                const Array4<Real>& sold_xmom = S_old[IntVar::xmom].array(mfi);
                const Array4<Real>& sold_ymom = S_old[IntVar::ymom].array(mfi);
                const Array4<Real>& sold_zmom = S_old[IntVar::zmom].array(mfi);

                ParallelFor(gtbx, gtby, gtbz,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                    scrh_xmom(i,j,k) = sold_xmom(i,j,k);
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                    scrh_ymom(i,j,k) = sold_ymom(i,j,k);
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                    scrh_zmom(i,j,k) = sold_zmom(i,j,k);
                });
            }
        }

        for (int nrk = 0; nrk < 3; nrk++)
        {
            // amrex::Print() << "Starting RK3: Step " << nrk+1 << std::endl;
            if (nrk == 0) { nsubsteps = 1;               dtau = timestep / 3.0; }
            if (nrk == 1) { nsubsteps = substep_ratio/2; dtau = sub_timestep;}
            if (nrk == 2) { nsubsteps = substep_ratio;   dtau = sub_timestep;}

            // step 1 starts with S_stage = S^n  and we always start substepping at the old time
            // step 2 starts with S_stage = S^*  and we always start substepping at the old time
            // step 3 starts with S_stage = S^** and we always start substepping at the old time

            // We fill the fast variables with the data from the start of the RK iteration, while
            // we fill the slow variables with the most recent RK stage data
            // Also, S_scratch will hold (rho theta) from the previous fast timestep
        #ifdef _OPENMP
        #pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
        #endif
            {
                for ( MFIter mfi(S_old[IntVar::cons],TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                    const Box gbx = mfi.tilebox().grow(S_old[IntVar::cons].nGrowVect());
                    const Box gtbx = mfi.nodaltilebox(0).grow(S_old[IntVar::xmom].nGrowVect());
                    const Box gtby = mfi.nodaltilebox(1).grow(S_old[IntVar::ymom].nGrowVect());
                    const Box gtbz = mfi.nodaltilebox(2).grow(S_old[IntVar::zmom].nGrowVect());
                    const Box gfbx = mfi.nodaltilebox(0).grow(S_old[IntVar::xflux].nGrowVect());
                    const Box gfby = mfi.nodaltilebox(1).grow(S_old[IntVar::yflux].nGrowVect());
                    const Box gfbz = mfi.nodaltilebox(2).grow(S_old[IntVar::zflux].nGrowVect());

                    Vector<Array4<Real> > sold_h(IntVar::NumVars);
                    Vector<Array4<Real> > snew_h(IntVar::NumVars);
                    Vector<Array4<Real> > ssum_h(IntVar::NumVars);
                    Vector<Array4<Real> > scrh_h(IntVar::NumVars);

                    for (int i = 0; i < IntVar::NumVars; ++i) {
                        ssum_h[i] = (*S_sum)[i].array(mfi);
                        sold_h[i] = S_old[i].array(mfi);
                        snew_h[i] = S_new[i].array(mfi);
                        scrh_h[i] = (*S_scratch)[i].array(mfi);
                    }

                    Gpu::AsyncVector<Array4<Real> > sold_d(IntVar::NumVars);
                    Gpu::AsyncVector<Array4<Real> > snew_d(IntVar::NumVars);
                    Gpu::AsyncVector<Array4<Real> > ssum_d(IntVar::NumVars);
                    Gpu::AsyncVector<Array4<Real> > scrh_d(IntVar::NumVars);

                    Gpu::copy(Gpu::hostToDevice, sold_h.begin(), sold_h.end(), sold_d.begin());
                    Gpu::copy(Gpu::hostToDevice, snew_h.begin(), snew_h.end(), snew_d.begin());
                    Gpu::copy(Gpu::hostToDevice, ssum_h.begin(), ssum_h.end(), ssum_d.begin());
                    Gpu::copy(Gpu::hostToDevice, scrh_h.begin(), scrh_h.end(), scrh_d.begin());

                    Array4<Real>* sold = sold_d.dataPtr();
                    Array4<Real>* snew = snew_d.dataPtr();
                    Array4<Real>* ssum = ssum_d.dataPtr();
                    Array4<Real>* scrh = scrh_d.dataPtr();

                    ParallelFor(gbx, Cons::NumVars,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
                        if (n >= scomp_fast[IntVar::cons] && n < scomp_fast[IntVar::cons] + ncomp_fast[IntVar::cons]) {
                            ssum[IntVar::cons](i,j,k,n) = sold[IntVar::cons](i,j,k,n);
                        }
                        if (n >= scomp_slow[IntVar::cons] && n < scomp_slow[IntVar::cons] + ncomp_slow[IntVar::cons]) {
                            ssum[IntVar::cons](i,j,k,n) = snew[IntVar::cons](i,j,k,n);
                        }
                        if (n >= scomp_rth[IntVar::cons] && n < scomp_rth[IntVar::cons] + ncomp_rth[IntVar::cons]) {
                            scrh[IntVar::cons](i,j,k,n) = sold[IntVar::cons](i,j,k,n);
                        }
                    });

                    ParallelFor(gtbx, gtby, gtbz,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                        ssum[IntVar::xmom](i,j,k) = sold[IntVar::xmom](i,j,k);
                    },
                    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                        ssum[IntVar::ymom](i,j,k) = sold[IntVar::ymom](i,j,k);
                    },
                    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                        ssum[IntVar::zmom](i,j,k) = sold[IntVar::zmom](i,j,k);
                    });

                    ParallelFor(
                    gfbx, Cons::NumVars,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept {
                        if (n >= scomp_fast[IntVar::xflux] && n < scomp_fast[IntVar::xflux] + ncomp_fast[IntVar::xflux]) {
                            ssum[IntVar::xflux](i,j,k,n) = sold[IntVar::xflux](i,j,k,n);
                        }
                        if (n >= scomp_slow[IntVar::xflux] && n < scomp_slow[IntVar::xflux] + ncomp_slow[IntVar::xflux]) {
                            ssum[IntVar::xflux](i,j,k,n) = snew[IntVar::xflux](i,j,k,n);
                        }
                    },
                    gfby, Cons::NumVars, 
                    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept {
                        if (n >= scomp_fast[IntVar::yflux] && n < scomp_fast[IntVar::yflux] + ncomp_fast[IntVar::yflux]) {
                            ssum[IntVar::yflux](i,j,k,n) = sold[IntVar::yflux](i,j,k,n);
                        }
                        if (n >= scomp_slow[IntVar::yflux] && n < scomp_slow[IntVar::yflux] + ncomp_slow[IntVar::yflux]) {
                            ssum[IntVar::yflux](i,j,k,n) = snew[IntVar::yflux](i,j,k,n);
                        }
                    },
                    gfbz, Cons::NumVars,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept {
                        if (n >= scomp_fast[IntVar::zflux] && n < scomp_fast[IntVar::zflux] + ncomp_fast[IntVar::zflux]) {
                            ssum[IntVar::zflux](i,j,k,n) = sold[IntVar::zflux](i,j,k,n);
                        }
                        if (n >= scomp_slow[IntVar::zflux] && n < scomp_slow[IntVar::zflux] + ncomp_slow[IntVar::zflux]) {
                            ssum[IntVar::zflux](i,j,k,n) = snew[IntVar::zflux](i,j,k,n);
                        }
                    });
                }
            }

            // S_scratch also holds the average momenta over the fast iterations --
            //    to be used to update the slow variables -- we will initialize with
            //    the momenta used in the first call to the slow_rhs, then update
            //    inside fast_rhs, then use these values in the later call to slow_rhs

            // Evaluate F_slow(S_stage) only for the fast variables
            slow_rhs(*F_slow, S_new, *S_scratch, time, RHSVar::fast);

            amrex::Real inv_fac = 1.0 / static_cast<amrex::Real>(nsubsteps);

            for (int ks = 0; ks < nsubsteps; ++ks)
            {
                // Evaluate F_fast(S_pert, S_old) on the fast variables
                // S_sum is used in the fast RHS to internally define S_pert = (S_sum - S_stage)
                fast_rhs(*F_pert, *F_slow, S_new, *S_sum, *S_scratch, dtau, inv_fac);

                // Update S_sum = S_pert + S_stage only for the fast variables
            #ifdef _OPENMP
            #pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
            #endif
                {
                    for ( MFIter mfi(S_old[IntVar::cons],TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                        const Box bx = mfi.tilebox();
                        const Box tbx = mfi.nodaltilebox(0);
                        const Box tby = mfi.nodaltilebox(1);
                        const Box tbz = mfi.nodaltilebox(2);

                        Vector<Array4<Real> > ssum_h(IntVar::NumVars);
                        Vector<Array4<Real> > fslow_h(IntVar::NumVars);
                        Vector<Array4<Real> > fpert_h(IntVar::NumVars);

                        for (int i = 0; i < IntVar::NumVars; ++i) {
                            ssum_h[i]  = (*S_sum)[i].array(mfi);
                            fslow_h[i] = (*F_slow)[i].array(mfi);
                            fpert_h[i] = (*F_pert)[i].array(mfi);
                        }

                        Gpu::AsyncVector<Array4<Real> > ssum_d(IntVar::NumVars);
                        Gpu::AsyncVector<Array4<Real> > fslow_d(IntVar::NumVars);
                        Gpu::AsyncVector<Array4<Real> > fpert_d(IntVar::NumVars);

                        Gpu::copy(Gpu::hostToDevice, ssum_h.begin(), ssum_h.end(), ssum_d.begin());
                        Gpu::copy(Gpu::hostToDevice, fslow_h.begin(), fslow_h.end(), fslow_d.begin());
                        Gpu::copy(Gpu::hostToDevice, fpert_h.begin(), fpert_h.end(), fpert_d.begin());

                        Array4<Real>* ssum = ssum_d.dataPtr();
                        Array4<Real>* fslow = fslow_d.dataPtr();
                        Array4<Real>* fpert = fpert_d.dataPtr();

                        ParallelFor(bx, ncomp_fast[IntVar::cons],
                        [=] AMREX_GPU_DEVICE (int i, int j, int k, int nn) {
                            const int n = scomp_fast[IntVar::cons] + nn;
                            ssum[IntVar::cons](i,j,k,n) += dtau * fslow[IntVar::cons](i,j,k,n);
                            ssum[IntVar::cons](i,j,k,n) += dtau * fpert[IntVar::cons](i,j,k,n);
                        });

                        ParallelFor(tbx, tby, tbz,
                        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                            ssum[IntVar::xmom](i,j,k) += dtau * fslow[IntVar::xmom](i,j,k);
                            ssum[IntVar::xmom](i,j,k) += dtau * fpert[IntVar::xmom](i,j,k);
                        },
                        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                            ssum[IntVar::ymom](i,j,k) += dtau * fslow[IntVar::ymom](i,j,k);
                            ssum[IntVar::ymom](i,j,k) += dtau * fpert[IntVar::ymom](i,j,k);
                        },
                        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                            ssum[IntVar::zmom](i,j,k) += dtau * fslow[IntVar::zmom](i,j,k);
                            ssum[IntVar::zmom](i,j,k) += dtau * fpert[IntVar::zmom](i,j,k);
                        });

                        ParallelFor(
                        tbx, ncomp_fast[IntVar::xflux],
                        [=] AMREX_GPU_DEVICE (int i, int j, int k, int nn) noexcept {
                            const int n = scomp_fast[IntVar::xflux] + nn;
                            ssum[IntVar::xflux](i,j,k,n) += dtau * fslow[IntVar::xflux](i,j,k,n);
                            ssum[IntVar::xflux](i,j,k,n) += dtau * fpert[IntVar::xflux](i,j,k,n);
                        },
                        tby, ncomp_fast[IntVar::yflux],
                        [=] AMREX_GPU_DEVICE (int i, int j, int k, int nn) noexcept {
                            const int n = scomp_fast[IntVar::yflux] + nn;
                            ssum[IntVar::yflux](i,j,k,n) += dtau * fslow[IntVar::yflux](i,j,k,n);
                            ssum[IntVar::yflux](i,j,k,n) += dtau * fpert[IntVar::yflux](i,j,k,n);
                        },
                        tbz, ncomp_fast[IntVar::zflux],
                        [=] AMREX_GPU_DEVICE (int i, int j, int k, int nn) noexcept {
                            const int n = scomp_fast[IntVar::zflux] + nn;
                            ssum[IntVar::zflux](i,j,k,n) += dtau * fslow[IntVar::zflux](i,j,k,n);
                            ssum[IntVar::zflux](i,j,k,n) += dtau * fpert[IntVar::zflux](i,j,k,n);
                        });
                    }
                }

                // Call the post-substep hook for S_sum at t = time + (ks+1) * sub_dt
                // Only call the post_substep on the fast variables
                int ngrow = 1;
                post_substep(*S_sum    , time + (ks+1) * dtau, ngrow, ngrow);

                // Copy only (rho theta), including ghost cells
                amrex::MultiFab::Copy((*S_scratch)[IntVar::cons], (*S_sum)[IntVar::cons], Cons::RhoTheta, Cons::RhoTheta, 1,
                                      (*S_sum)[IntVar::cons].nGrowVect());
            }

            // Evaluate F_slow(S_stage) only for the slow variables
            // Note that we are using the current stage versions (in S_new) of the slow variables
            //      (because we didn't update the slow variables in the substepping)
            //       but we are using the "new" versions (in S_sum) of the velocities
            //      (because we did    update the fast variables in the substepping)
            slow_rhs(*F_slow, *S_sum, *S_scratch, time, RHSVar::slow);

            // Apply F_slow only for the slow variables
            // Then copy into S_new which holds the stage data
        #ifdef _OPENMP
        #pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
        #endif
            {
                for ( MFIter mfi(S_old[IntVar::cons],TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                    const Box bx = mfi.tilebox();
                    const Box tbx = mfi.nodaltilebox(0);
                    const Box tby = mfi.nodaltilebox(1);
                    const Box tbz = mfi.nodaltilebox(2);
                    const Box gbx = mfi.tilebox().grow(S_old[IntVar::cons].nGrowVect());
                    const Box gtbx = mfi.nodaltilebox(0).grow(S_old[IntVar::xmom].nGrowVect());
                    const Box gtby = mfi.nodaltilebox(1).grow(S_old[IntVar::ymom].nGrowVect());
                    const Box gtbz = mfi.nodaltilebox(2).grow(S_old[IntVar::zmom].nGrowVect());
                    const Box gfbx = mfi.nodaltilebox(0).grow(S_old[IntVar::xflux].nGrowVect());
                    const Box gfby = mfi.nodaltilebox(1).grow(S_old[IntVar::yflux].nGrowVect());
                    const Box gfbz = mfi.nodaltilebox(2).grow(S_old[IntVar::zflux].nGrowVect());

                    Vector<Array4<Real> > sold_h(IntVar::NumVars);
                    Vector<Array4<Real> > snew_h(IntVar::NumVars);
                    Vector<Array4<Real> > ssum_h(IntVar::NumVars);
                    Vector<Array4<Real> > fslow_h(IntVar::NumVars);

                    for (int i = 0; i < IntVar::NumVars; ++i) {
                        ssum_h[i]  = (*S_sum)[i].array(mfi);
                        fslow_h[i] = (*F_slow)[i].array(mfi);
                        sold_h[i] = S_old[i].array(mfi);
                        snew_h[i] = S_new[i].array(mfi);
                    }

                    Gpu::AsyncVector<Array4<Real> > sold_d(IntVar::NumVars);
                    Gpu::AsyncVector<Array4<Real> > snew_d(IntVar::NumVars);
                    Gpu::AsyncVector<Array4<Real> > ssum_d(IntVar::NumVars);
                    Gpu::AsyncVector<Array4<Real> > fslow_d(IntVar::NumVars);

                    Gpu::copy(Gpu::hostToDevice, sold_h.begin(), sold_h.end(), sold_d.begin());
                    Gpu::copy(Gpu::hostToDevice, snew_h.begin(), snew_h.end(), snew_d.begin());
                    Gpu::copy(Gpu::hostToDevice, ssum_h.begin(), ssum_h.end(), ssum_d.begin());
                    Gpu::copy(Gpu::hostToDevice, fslow_h.begin(), fslow_h.end(), fslow_d.begin());

                    Array4<Real>* sold = sold_d.dataPtr();
                    Array4<Real>* snew = snew_d.dataPtr();
                    Array4<Real>* ssum = ssum_d.dataPtr();
                    Array4<Real>* fslow = fslow_d.dataPtr();

                    ParallelFor(gbx, ncomp_slow[IntVar::cons],
                    [=] AMREX_GPU_DEVICE (int i, int j, int k, int nn) noexcept {
                        const int n = scomp_slow[IntVar::cons] + nn;
                        ssum[IntVar::cons](i,j,k,n) = sold[IntVar::cons](i,j,k,n);
                        if (bx.contains(i,j,k))
                            ssum[IntVar::cons](i,j,k,n) += nsubsteps * dtau * fslow[IntVar::cons](i,j,k,n);
                    });

                    ParallelFor(gbx, ncomp_all[IntVar::cons],
                    [=] AMREX_GPU_DEVICE (int i, int j, int k, int nn) noexcept {
                        const int n = scomp_all[IntVar::cons] + nn;
                        snew[IntVar::cons](i,j,k,n) = ssum[IntVar::cons](i,j,k,n);
                    });

                    ParallelFor(gtbx, gtby, gtbz,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                        snew[IntVar::xmom](i,j,k) = ssum[IntVar::xmom](i,j,k);
                    },
                    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                        snew[IntVar::ymom](i,j,k) = ssum[IntVar::ymom](i,j,k);
                    },
                    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                        snew[IntVar::zmom](i,j,k) = ssum[IntVar::zmom](i,j,k);
                    });

                    ParallelFor(
                    gfbx, ncomp_slow[IntVar::xflux],
                    [=] AMREX_GPU_DEVICE (int i, int j, int k, int nn) noexcept {
                        const int n = scomp_slow[IntVar::xflux] + nn;
                        ssum[IntVar::xflux](i,j,k,n) = sold[IntVar::xflux](i,j,k,n);
                        if (tbx.contains(i,j,k))
                            ssum[IntVar::xflux](i,j,k,n) += nsubsteps * dtau * fslow[IntVar::xflux](i,j,k,n);
                    },
                    gfby, ncomp_slow[IntVar::yflux],
                    [=] AMREX_GPU_DEVICE (int i, int j, int k, int nn) noexcept {
                        const int n = scomp_slow[IntVar::yflux] + nn;
                        ssum[IntVar::yflux](i,j,k,n) = sold[IntVar::yflux](i,j,k,n);
                        if (tby.contains(i,j,k))
                            ssum[IntVar::yflux](i,j,k,n) += nsubsteps * dtau * fslow[IntVar::yflux](i,j,k,n);
                    },
                    gfbz, ncomp_slow[IntVar::zflux],
                    [=] AMREX_GPU_DEVICE (int i, int j, int k, int nn) noexcept {
                        const int n = scomp_slow[IntVar::zflux] + nn;
                        ssum[IntVar::zflux](i,j,k,n) = sold[IntVar::zflux](i,j,k,n);
                        if (tbz.contains(i,j,k))
                            ssum[IntVar::zflux](i,j,k,n) += nsubsteps * dtau * fslow[IntVar::zflux](i,j,k,n);
                    });

                    ParallelFor(
                    gfbx, ncomp_all[IntVar::xflux],
                    [=] AMREX_GPU_DEVICE (int i, int j, int k, int nn) noexcept {
                        const int n = scomp_all[IntVar::xflux] + nn;
                        snew[IntVar::xflux](i,j,k,n) = ssum[IntVar::xflux](i,j,k,n);
                    },
                    gfby, ncomp_all[IntVar::yflux],
                    [=] AMREX_GPU_DEVICE (int i, int j, int k, int nn) noexcept {
                        const int n = scomp_all[IntVar::yflux] + nn;
                        snew[IntVar::yflux](i,j,k,n) = ssum[IntVar::yflux](i,j,k,n);
                    },
                    gfbz, ncomp_all[IntVar::zflux],
                    [=] AMREX_GPU_DEVICE (int i, int j, int k, int nn) noexcept {
                        const int n = scomp_all[IntVar::zflux] + nn;
                        snew[IntVar::zflux](i,j,k,n) = ssum[IntVar::zflux](i,j,k,n);
                    });
                }
            }

            // Call the post-update hook for S_new after all the fine steps completed
            // This will update S_prim that is used in the slow RHS
            if (nrk < 2)
               post_update(S_new, time + nsubsteps*dtau,S_new[IntVar::cons].nGrow(), S_new[IntVar::xmom].nGrow());
        }

        // Return timestep
        return timestep;
    }

    void time_interpolate (const T& /* S_new */, const T& /*S_old*/, amrex::Real /*timestep_fraction*/, T& /*data*/) {}

    void map_data (std::function<void(T&)> Map)
    {
        for (auto& F : T_store) {
            Map(*F);
        }
    }

};

#endif
