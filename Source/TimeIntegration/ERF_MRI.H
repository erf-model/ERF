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
    std::function<void(T&, const T&,     const amrex::Real, const amrex::Real     )> rhs;
    std::function<void(T&, T&, T&, const amrex::Real, const amrex::Real     )> slow_rhs_pre;
    std::function<void(T&, T&, T&, T&, T&, const amrex::Real, const amrex::Real     )> slow_rhs_post;
    std::function<void(int, T&, const T&, T&, T&, T&, const amrex::Real, const amrex::Real, const amrex::Real)> fast_rhs;

   /**
    * \brief Integrator timestep size (Real)
    */
    amrex::Real timestep;

   /**
    * \brief The ratio of slow timestep size / fast timestep size (int)
    */
    int slow_fast_timestep_ratio = 0;

   /**
    * \brief Should we not do acoustic substepping
    */
    int no_substepping;

   /**
    * \brief The post_update function is called by the integrator on state data before using it to evaluate a right-hand side.
    */
    std::function<void (T&, amrex::Real, amrex::Real, int, int)> post_update;
    std::function<void (T&, T&, T&, amrex::Real, amrex::Real)> no_substep;


    amrex::Vector<std::unique_ptr<T> > T_store;
    T* S_sum;
    T* S_scratch;
    T* F_slow;

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

    void setNoSubstepping(int _no_substepping)
    {
        no_substepping = _no_substepping;
    }

    void set_rhs (std::function<void(T&, const T&, const amrex::Real, const amrex::Real, int)> F)
    {
        rhs = F;
    }

    void set_slow_rhs_pre (std::function<void(T&, T&, T&, const amrex::Real, const amrex::Real)> F)
    {
        slow_rhs_pre = F;
    }
    void set_slow_rhs_post (std::function<void(T&, T&, T&, T&, T&, const amrex::Real, const amrex::Real)> F)
    {
        slow_rhs_post = F;
    }

    void set_fast_rhs (std::function<void(int, T&, const T&, T&, T&, T&, const amrex::Real, const amrex::Real, const amrex::Real)> F)
    {
        fast_rhs = F;
    }

    void set_slow_fast_timestep_ratio (const int timestep_ratio = 1)
    {
        slow_fast_timestep_ratio = timestep_ratio;
    }

    int get_slow_fast_timestep_ratio ()
    {
        return slow_fast_timestep_ratio;
    }

    void set_post_update (std::function<void (T&, amrex::Real, amrex::Real, int, int)> F)
    {
        post_update = F;
    }

    void set_no_substep (std::function<void (T&, T&, T&, amrex::Real, amrex::Real)> F)
    {
        no_substep = F;
    }

    std::function<void(T&, const T&, const amrex::Real, int)> get_rhs ()
    {
        return rhs;
    }

    amrex::Real advance (T& S_old, T& S_new, amrex::Real time, const amrex::Real time_step)
    {
    BL_PROFILE_REGION("MRI_advance");
        using namespace amrex;

        // *******************************************************************************
        // version == 0: we only update the fast variables every fast timestep, then update
        //               the slow variables after the acoustic sub-stepping.  This has
        //               two calls to slow_rhs so that we can update the slow variables
        //               with the velocity field after the acoustic substepping using
        //               the time-averaged velocity from the substepping
        // version == 1: we don't do any acoustic subcyling so we only make one call per RK
        //               stage to slow_rhs
        // *******************************************************************************
        int version = no_substepping ? 1 : 0;

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

        int n_data = IntVar::NumVars;

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

                Vector<Array4<Real> > sold_h(n_data);
                Vector<Array4<Real> > snew_h(n_data);

                for (int i = 0; i < n_data; ++i) {
                    sold_h[i] = S_old[i].array(mfi);
                    snew_h[i] = S_new[i].array(mfi);
                }

                Gpu::AsyncVector<Array4<Real> > sold_d(n_data);
                Gpu::AsyncVector<Array4<Real> > snew_d(n_data);

                Gpu::copy(Gpu::hostToDevice, sold_h.begin(), sold_h.end(), sold_d.begin());
                Gpu::copy(Gpu::hostToDevice, snew_h.begin(), snew_h.end(), snew_d.begin());

                Array4<Real>* sold = sold_d.dataPtr();
                Array4<Real>* snew = snew_d.dataPtr();

                ParallelFor(gbx, static_cast<int>(Cons::NumVars),
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
            }
        }

        post_update(S_new, time, time_step, S_new[IntVar::cons].nGrow(), S_new[IntVar::xmom].nGrow());

        // Timestep taken by the fast integrator
        amrex::Real dtau;

        // How many timesteps taken by the fast integrator
        int nsubsteps;

        int nav = Cons::NumVars;
        const amrex::GpuArray<int, IntVar::NumVars> scomp_all = {0,0,0,0};
        const amrex::GpuArray<int, IntVar::NumVars> ncomp_all = {nav,1,1,1};

        const amrex::GpuArray<int, IntVar::NumVars> scomp_fast = {0,0,0,0};
        const amrex::GpuArray<int, IntVar::NumVars> ncomp_fast = {2,1,1,1};

        int nsv = Cons::NumVars-2;
        const amrex::GpuArray<int, IntVar::NumVars> scomp_slow = {  2,0,0,0};
        const amrex::GpuArray<int, IntVar::NumVars> ncomp_slow = {nsv,0,0,0};

        const amrex::GpuArray<int, IntVar::NumVars> scomp_rth  = {1,0,0,0};
        const amrex::GpuArray<int, IntVar::NumVars> ncomp_rth  = {1,0,0,0};

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
            amrex::Real time_stage;
            if (nrk == 0) { nsubsteps = 1;               dtau = timestep / 3.0; time_stage = time + timestep / 3.0;}
            if (nrk == 1) { nsubsteps = substep_ratio/2; dtau = sub_timestep  ; time_stage = time + timestep / 2.0;}
            if (nrk == 2) { nsubsteps = substep_ratio;   dtau = sub_timestep  ; time_stage = time + timestep      ;}

            // step 1 starts with S_stage = S^n  and we always start substepping at the old time
            // step 2 starts with S_stage = S^*  and we always start substepping at the old time
            // step 3 starts with S_stage = S^** and we always start substepping at the old time

#if 0
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

                    Vector<Array4<Real> > sold_h(n_data);
                    Vector<Array4<Real> > snew_h(n_data);
                    Vector<Array4<Real> > ssum_h(n_data);
                    Vector<Array4<Real> > scrh_h(n_data);

                    for (int i = 0; i < n_data; ++i) {
                        ssum_h[i] = (*S_sum)[i].array(mfi);
                        sold_h[i] = S_old[i].array(mfi);
                        snew_h[i] = S_new[i].array(mfi);
                        scrh_h[i] = (*S_scratch)[i].array(mfi);
                    }

                    Gpu::AsyncVector<Array4<Real> > sold_d(n_data);
                    Gpu::AsyncVector<Array4<Real> > snew_d(n_data);
                    Gpu::AsyncVector<Array4<Real> > ssum_d(n_data);
                    Gpu::AsyncVector<Array4<Real> > scrh_d(n_data);

                    Gpu::copy(Gpu::hostToDevice, sold_h.begin(), sold_h.end(), sold_d.begin());
                    Gpu::copy(Gpu::hostToDevice, snew_h.begin(), snew_h.end(), snew_d.begin());
                    Gpu::copy(Gpu::hostToDevice, ssum_h.begin(), ssum_h.end(), ssum_d.begin());
                    Gpu::copy(Gpu::hostToDevice, scrh_h.begin(), scrh_h.end(), scrh_d.begin());

                    Array4<Real>* sold = sold_d.dataPtr();
                    Array4<Real>* snew = snew_d.dataPtr();
                    Array4<Real>* ssum = ssum_d.dataPtr();
                    Array4<Real>* scrh = scrh_d.dataPtr();

                    if (version == 0) {
                        ParallelFor(gbx, static_cast<int>(Cons::NumVars),
                        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
//                          if (n >= scomp_slow[IntVar::cons] && n < scomp_slow[IntVar::cons] + ncomp_slow[IntVar::cons]) {
//                              ssum[IntVar::cons](i,j,k,n) = snew[IntVar::cons](i,j,k,n);
//                          }
                            if (n >= scomp_rth[IntVar::cons] && n < scomp_rth[IntVar::cons] + ncomp_rth[IntVar::cons]) {
                                scrh[IntVar::cons](i,j,k,n) = sold[IntVar::cons](i,j,k,n);
                            }
                        });
                    } // with substepping only
                }
            }
#endif

            // S_scratch also holds the average momenta over the fast iterations --
            //    to be used to update the slow variables -- we will initialize with
            //    the momenta used in the first call to the slow_rhs, then update
            //    inside fast_rhs, then use these values in the later call to slow_rhs

            // ****************************************************
            // if (version == 0) evaluate F_slow(S_stage) only for the fast variables
            // if (version == 1) evaluate F_slow(S_stage) only for all the  variables
            // ****************************************************
            slow_rhs_pre(*F_slow, S_new, *S_scratch, time, time_stage);

            amrex::Real inv_fac = 1.0 / static_cast<amrex::Real>(nsubsteps);

            // ****************************************************
            // Acoustic substepping
            // ****************************************************
            if (version == 0)
            {
                // *******************************************************************************
                // Update the fast variables
                // *******************************************************************************
                for (int ks = 0; ks < nsubsteps; ++ks)
                {
                    fast_rhs(ks, *F_slow, S_old, S_new, *S_sum, *S_scratch, dtau, inv_fac, time + (ks+1) * dtau);

                } // ks

            } else {
                no_substep(*S_sum, S_old, *F_slow, time + nsubsteps*dtau, nsubsteps*dtau);
            }

            // ****************************************************
            // Evaluate F_slow(S_stage) only for the slow variables
            // Note that we are using the current stage versions (in S_new) of the slow variables
            //      (because we didn't update the slow variables in the substepping)
            //       but we are using the "new" versions (in S_sum) of the velocities
            //      (because we did    update the fast variables in the substepping)
            // ****************************************************
            slow_rhs_post(*F_slow, S_old, S_new, *S_sum, *S_scratch, time, time_stage);

            // Call the post-update hook for S_new after all the fast steps completed
            // This will update S_prim that is used in the slow RHS
            post_update(S_new, time + nsubsteps*dtau, dtau, S_new[IntVar::cons].nGrow(), S_new[IntVar::xmom].nGrow());
        } // nrk

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
