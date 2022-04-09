#ifndef ERF_MRI_H
#define ERF_MRI_H
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_ParmParse.H>
#include <AMReX_IntegratorBase.H>
#include <TimeIntegration.H>
#include <functional>

using namespace amrex;

template<class T>
class MRISplitIntegrator : public IntegratorBase<T>
{
private:
    typedef IntegratorBase<T> BaseT;

    std::function<void (T&, amrex::Real)> post_substep;
    std::function<void(T&, T&, T&, const T&, const T&, const amrex::Real, const amrex::Real)> implicit_fast_rhs;

    amrex::Vector<std::unique_ptr<T> > T_store;
    T* S_sum;
    T* S_sum_old;
    T* F_slow;
    T* F_pert;

    void initialize_data (const T& S_data)
    {
        const bool include_ghost = true;
        IntegratorOps<T>::CreateLike(T_store, S_data, include_ghost);
        S_sum = T_store[0].get();
        IntegratorOps<T>::CreateLike(T_store, S_data, include_ghost);
        S_sum_old = T_store[1].get();
        IntegratorOps<T>::CreateLike(T_store, S_data, include_ghost);
        F_slow = T_store[2].get();
        IntegratorOps<T>::CreateLike(T_store, S_data, include_ghost);
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

    void set_post_substep (std::function<void (T&, amrex::Real)> F)
    {
        post_substep = F;
    }

    std::function<void (T&, amrex::Real)> get_post_substep ()
    {
        return post_substep;
    }

    void set_implicit_fast_rhs (std::function<void(T&, T&, T&, const T&, const T&, const amrex::Real, const amrex::Real)> F)
    {
        implicit_fast_rhs = F;
    }

    std::function<void(T&, T&, T&, const T&, const amrex::Real, const amrex::Real)> get_implicit_fast_rhs ()
    {
        return implicit_fast_rhs;
    }

    amrex::Real advance (T& S_old, T& S_new, amrex::Real time, const amrex::Real time_step)
    {
        BaseT::timestep = time_step;

        const int substep_ratio = BaseT::get_slow_fast_timestep_ratio();

        AMREX_ALWAYS_ASSERT(substep_ratio > 1 && substep_ratio % 2 == 0);

        const Real sub_timestep = BaseT::timestep / substep_ratio;

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
        /* RK3 Integration with Acoustuc Sub-stepping */
        /**********************************************/

        // Start with S_new (aka S_stage) holding S_old
        IntegratorOps<T>::Copy(S_new, S_old);
        BaseT::post_update(S_new, time);

        // Timestep taken by the fast integrator
        amrex::Real dtau;

        // How many timesteps taken by the fast integrator
        int nsubsteps;

        for (int nrk = 0; nrk < 3; nrk++)
        {
            amrex::Print() << "Starting RK3: Step " << nrk+1 << std::endl;
            if (nrk == 0) { nsubsteps = 1;               dtau = BaseT::timestep / 3.0_rt; }
            if (nrk == 1) { nsubsteps = substep_ratio/2; dtau = sub_timestep;}
            if (nrk == 2) { nsubsteps = substep_ratio;   dtau = sub_timestep;}

            // step 1 starts with S_stage = S^n  and we always start substepping at the old time
            // step 2 starts with S_stage = S^*  and we always start substepping at the old time
            // step 3 starts with S_stage = S^** and we always start substepping at the old time
            IntegratorOps<T>::Copy(*S_sum    , S_old);
            IntegratorOps<T>::Copy(*S_sum_old, S_old);

            // Evaluate F_slow(S_stage)
            BaseT::rhs(*F_slow, S_new, time);

            for (int k = 0; k < nsubsteps; ++k)
            {
                // Evaluate F_pert(S_pert, S_old)
                // S_sum is used in the fast RHS to internally define S_pert = (S_sum - S_stage)
                implicit_fast_rhs(*F_pert, *F_slow, S_new, *S_sum, *S_sum_old, time + k*dtau, dtau);

                // Update S_sum = S_pert + S_stage
                IntegratorOps<T>::Saxpy(*S_sum, dtau, *F_slow);
                IntegratorOps<T>::Saxpy(*S_sum, dtau, *F_pert);

                IntegratorOps<T>::Copy(*S_sum_old, *S_sum);

                // Call the post-substep hook for S_sum at t = time + (k+1) * sub_dt
                post_substep(*S_sum    , time + (k+1) * dtau);
                post_substep(*S_sum_old, time + (k+1) * dtau);
            }

            // Define the solution at the next stage
            IntegratorOps<T>::Copy(S_new, *S_sum);

            // Call the post-update hook for S_new after all the fine steps completed
            BaseT::post_update(S_new, time + nsubsteps*dtau);
        }

        // Return timestep
        return BaseT::timestep;
    }

    void time_interpolate (const T& /* S_new */, const T& S_old, amrex::Real timestep_fraction, T& data) {}

    void map_data (std::function<void(T&)> Map)
    {
        for (auto& F : T_store) {
            Map(*F);
        }
    }

};

#endif
