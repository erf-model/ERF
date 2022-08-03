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
        amrex::IntegratorOps<T>::Copy(S_new, S_old);
        post_update(S_new, time, S_new[IntVar::cons].nGrow(), S_new[IntVar::xmom].nGrow());

        // Timestep taken by the fast integrator
        amrex::Real dtau;

        // How many timesteps taken by the fast integrator
        int nsubsteps;

        int nav = Cons::NumVars;
        const amrex::Vector<int> scomp_all = {0,0,0,0,0,0,0};
        const amrex::Vector<int> ncomp_all = {nav,1,1,1,nav,nav,nav};

        const amrex::Vector<int> scomp_fast = {0,0,0,0,0,0,0};
        const amrex::Vector<int> ncomp_fast = {2,1,1,1,2,2,2};

        int nsv = Cons::NumVars-2;
        const amrex::Vector<int> scomp_slow = {  2,0,0,0,  2,  2,  2};
        const amrex::Vector<int> ncomp_slow = {nsv,0,0,0,nsv,nsv,nsv};

        // We copy (rho theta) and the velocities here -- the velocities
        //    will be over-written in slow_rhs on all valid faces but we
        //    use this copy to fill in the ghost locations which will
        //    be needed for metric terms
        const amrex::Vector<int> scomp_vels = {0,0,0,0,0,0,0};
        const amrex::Vector<int> ncomp_vels = {0,1,1,1,0,0,0};
        amrex::IntegratorOps<T>::Copy(*S_scratch, S_old, scomp_vels, ncomp_vels);

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
            amrex::IntegratorOps<T>::Copy(*S_sum    , S_old, scomp_fast, ncomp_fast);
            amrex::IntegratorOps<T>::Copy(*S_sum    , S_new, scomp_slow, ncomp_slow);

            // S_scratch will hold (rho theta) from the previous fast timestep
            const amrex::Vector<int> scomp_rth  = {1,0,0,0,0,0,0};
            const amrex::Vector<int> ncomp_rth  = {1,0,0,0,0,0,0};
            amrex::IntegratorOps<T>::Copy(*S_scratch, S_old, scomp_rth, ncomp_rth);

            // S_scratch also holds the average momenta over the fast iterations --
            //    to be used to update the slow variables -- we will initialize with
            //    the momenta used in the first call to the slow_rhs, then update
            //    inside fast_rhs, then use these values in the later call to slow_rhs

            // Evaluate F_slow(S_stage) only for the fast variables
            slow_rhs(*F_slow, S_new, *S_scratch, time, RHSVar::fast);

            amrex::Real inv_fac = 1.0 / static_cast<amrex::Real>(nsubsteps);

            for (int k = 0; k < nsubsteps; ++k)
            {
                // Evaluate F_fast(S_pert, S_old) on the fast variables
                // S_sum is used in the fast RHS to internally define S_pert = (S_sum - S_stage)
                fast_rhs(*F_pert, *F_slow, S_new, *S_sum, *S_scratch, dtau, inv_fac);

                // Update S_sum = S_pert + S_stage only for the fast variables
                amrex::IntegratorOps<T>::Saxpy(*S_sum, dtau, *F_slow, scomp_fast, ncomp_fast);
                amrex::IntegratorOps<T>::Saxpy(*S_sum, dtau, *F_pert, scomp_fast, ncomp_fast);

                // Call the post-substep hook for S_sum at t = time + (k+1) * sub_dt
                // Only call the post_substep on the fast variables
                int ngrow = 1;
                post_substep(*S_sum    , time + (k+1) * dtau, ngrow, ngrow);

                // Copy only (rho theta), including ghost cells
                amrex::IntegratorOps<T>::Copy(*S_scratch, *S_sum, scomp_rth, ncomp_rth, true);
            }
#if 0
            // Update S_sum = S_pert + S_stage only for the fast variables
            amrex::IntegratorOps<T>::Saxpy(*S_sum, time_step, *F_slow, scomp_fast, ncomp_fast);
#endif

            // Evaluate F_slow(S_stage) only for the slow variables
            // Note that we are using the current stage versions (in S_new) of the slow variables
            //      (because we didn't update the slow variables in the substepping)
            //       but we are using the "new" versions (in S_sum) of the velocities
            //      (because we did    update the fast variables in the substepping)
            slow_rhs(*F_slow, *S_sum, *S_scratch, time, RHSVar::slow);

            // Apply F_slow only for the slow variables
            amrex::IntegratorOps<T>::Copy(*S_sum    , S_old, scomp_slow, ncomp_slow);
            amrex::IntegratorOps<T>::Saxpy(*S_sum, nsubsteps*dtau, *F_slow, scomp_slow, ncomp_slow);

            // Copy into S_new which holds the stage data
            amrex::IntegratorOps<T>::Copy(S_new, *S_sum, scomp_all, ncomp_all);

            // Call the post-update hook for S_new after all the fine steps completed
            // This will update S_prim that is used in the slow RHS
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
