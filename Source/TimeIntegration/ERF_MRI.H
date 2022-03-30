#ifndef ERF_MRI_H
#define ERF_MRI_H
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_ParmParse.H>
#include <AMReX_IntegratorBase.H>
#include <functional>

using namespace amrex;

template<class T>
class MRISplitIntegrator : public IntegratorBase<T>
{
private:
    typedef IntegratorBase<T> BaseT;

    std::function<void (T&, amrex::Real)> post_substep;
    std::function<void(T&, T&, T&, const T&, const amrex::Real, const amrex::Real)> implicit_fast_rhs;
    bool use_implicit_fast_rhs;

    amrex::Vector<std::unique_ptr<T> > T_store;
    T* S_pert;
    T* S_sum;
    T* F_bar;
    T* F_pert;

    void initialize_data (const T& S_data)
    {
        const bool include_ghost = true;
        IntegratorOps<T>::CreateLike(T_store, S_data, include_ghost);
        S_pert = T_store[0].get();
        IntegratorOps<T>::CreateLike(T_store, S_data, include_ghost);
        S_sum = T_store[1].get();
        IntegratorOps<T>::CreateLike(T_store, S_data, include_ghost);
        F_bar = T_store[2].get();
        IntegratorOps<T>::CreateLike(T_store, S_data, include_ghost);
        F_pert = T_store[3].get();

    }

public:
    MRISplitIntegrator () : use_implicit_fast_rhs(false) {}

    MRISplitIntegrator (const T& S_data) : use_implicit_fast_rhs(false)
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

    void set_implicit_fast_rhs (std::function<void(T&, T&, T&, const T&, const amrex::Real, const amrex::Real)> F)
    {
        implicit_fast_rhs = F;
        use_implicit_fast_rhs = true;
    }

    std::function<void(T&, T&, T&, const T&, const amrex::Real, const amrex::Real)> get_implicit_fast_rhs ()
    {
        return implicit_fast_rhs;
    }

    amrex::Real advance (T& S_old, T& S_new, amrex::Real time, const amrex::Real time_step)
    {
        BaseT::timestep = time_step;

        const int substep_ratio = BaseT::get_slow_fast_timestep_ratio();
        AMREX_ALWAYS_ASSERT(substep_ratio >= 1 && substep_ratio % 2 == 0);
        const Real sub_timestep = BaseT::timestep / substep_ratio;

        // Assume before advance() that S_old is valid data at the current time ("time" argument)
        // And that if data is a MultiFab, both S_old and S_new contain ghost cells for evaluating a stencil based RHS
        // We need this from S_old. This is convenient for S_new to have so we can use it
        // as scratch space for stage values without creating a new scratch MultiFab with ghost cells.

        // NOTE: In the following, we use S_new to hold S*, S**, and finally, S^(n+1) at the new time
        // DEFINITIONS:
        // S_pert = S(t) - S_bar
        // S_sum = S_pert + S_bar = S(t)
        // F_bar = F(S_bar)
        // F_pert = G(S_pert, S_bar)

        /***********/
        /* STEP 1: */
        /***********/

        // S_pert = 0 since step 1 starts with S_bar = S_old
        IntegratorOps<T>::Copy(*S_pert, S_old);
        IntegratorOps<T>::Saxpy(*S_pert, -1.0_rt, S_old);

        // S_sum = S_pert + S_bar = 0 + S_old = S_old
        IntegratorOps<T>::Copy(*S_sum, S_old);

        // Evaluate F_bar(S_old)
        BaseT::rhs(*F_bar, S_old, time);

        // Evaluate F_pert(S_pert, S_old)
        if (use_implicit_fast_rhs)
            implicit_fast_rhs(*F_pert, *F_bar, S_old, *S_sum, time, BaseT::timestep / 3.0_rt);
        else
            BaseT::fast_rhs(*F_pert, S_old, *S_sum, time);

        // Calculate S* = S_old + (dt/3) * (Fbar(S_old) + Fpert(S_pert, S_old))
        IntegratorOps<T>::Copy(S_new, S_old);
        IntegratorOps<T>::Saxpy(S_new, BaseT::timestep / 3.0_rt, *F_bar);
        IntegratorOps<T>::Saxpy(S_new, BaseT::timestep / 3.0_rt, *F_pert);

        // Call the post-update hook for S_new=S* at t = time + dt/3
        BaseT::post_update(S_new, time + BaseT::timestep/3.0_rt);

        /***********/
        /* STEP 2: */
        /***********/

        // step 2 starts with S_bar = S* and we start substepping at the old time so S_pert = S_old - S*
        IntegratorOps<T>::Copy(*S_pert, S_old);
        IntegratorOps<T>::Saxpy(*S_pert, -1.0_rt, S_new);

        // S_sum = S_pert + S_bar = (S_old - S*) + S* = S_old
        IntegratorOps<T>::Copy(*S_sum, S_old);

        // Calculate F_bar = F(S_bar) = F(S*) at t = time + dt/3
        BaseT::rhs(*F_bar, S_new, time + BaseT::timestep / 3.0_rt);

        // Do forward Euler substepping for the fast RHS
        // from t = time to time + dt/2
        for (int k = 0; k < substep_ratio/2; ++k) {
            // S_pert^(k+1) = S_pert^(k) + dtau * (F_bar + F_pert)
            // first, evaluate F_pert(S_pert, S_bar)
            if (use_implicit_fast_rhs)
                implicit_fast_rhs(*F_pert, *F_bar, S_new, *S_sum, time + k * sub_timestep, sub_timestep);
            else
                BaseT::fast_rhs(*F_pert, S_new, *S_sum, time + k * sub_timestep);

            // now, get the new S_pert
            IntegratorOps<T>::Saxpy(*S_pert, sub_timestep, *F_bar);
            IntegratorOps<T>::Saxpy(*S_pert, sub_timestep, *F_pert);

            // next, recalculate S_sum = S_pert + S_bar
            IntegratorOps<T>::Copy(*S_sum, S_new);
            IntegratorOps<T>::Saxpy(*S_sum, 1.0_rt, *S_pert);

            // Call the post-substep hook for S_sum at t = time + (k+1) * sub_dt
            post_substep(*S_sum, time + (k+1) * sub_timestep);
        }

        // calculate S** = S* + S_pert
        IntegratorOps<T>::Saxpy(S_new, 1.0_rt, *S_pert);

        // Call the post-update hook for S_new=S** at t = time + dt/2
        BaseT::post_update(S_new, time + BaseT::timestep/2.0_rt);

        /***********/
        /* STEP 3: */
        /***********/

        // step 3 starts with S_bar = S** and we start substepping at the old time so S_pert = S_old - S**
        IntegratorOps<T>::Copy(*S_pert, S_old);
        IntegratorOps<T>::Saxpy(*S_pert, -1.0_rt, S_new);

        // S_sum = S_pert + S_bar = (S_old - S**) + S** = S_old
        IntegratorOps<T>::Copy(*S_sum, S_old);

        // Calculate F_bar = F(S_bar) = F(S**) at t = time + dt/2
        BaseT::rhs(*F_bar, S_new, time + BaseT::timestep / 2.0_rt);

        // Do forward Euler substepping for the fast RHS
        // from t = time to time + dt
        for (int k = 0; k < substep_ratio; ++k) {
            // S_pert^(k+1) = S_pert^(k) + dtau * (F_bar + F_pert)
            // first, evaluate F_pert(S_pert, S_bar)
            if (use_implicit_fast_rhs)
                implicit_fast_rhs(*F_pert, *F_bar, S_new, *S_sum, time + k * sub_timestep, sub_timestep);
            else
                BaseT::fast_rhs(*F_pert, S_new, *S_sum, time + k * sub_timestep);

            // now, get the new S_pert
            IntegratorOps<T>::Saxpy(*S_pert, sub_timestep, *F_bar);
            IntegratorOps<T>::Saxpy(*S_pert, sub_timestep, *F_pert);

            // next, recalculate S_sum = S_pert + S_bar
            IntegratorOps<T>::Copy(*S_sum, S_new);
            IntegratorOps<T>::Saxpy(*S_sum, 1.0_rt, *S_pert);

            // Call the post-substep hook for S_sum at t = time + (k+1) * sub_dt
            post_substep(*S_sum, time + (k+1) * sub_timestep);
        }

        // calculate S^(n+1) = S** + S_pert
        IntegratorOps<T>::Saxpy(S_new, 1.0_rt, *S_pert);

        // Call the post-update hook for S_new = S^(n+1) at t = time + dt
        BaseT::post_update(S_new, time + BaseT::timestep);

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
