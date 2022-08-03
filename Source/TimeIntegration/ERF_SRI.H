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
    std::function<void(T&, const T&, T&, const amrex::Real, int)> rhs;

   /**
    * \brief Integrator timestep size (Real)
    */
    amrex::Real timestep;

   /**
    * \brief The post_update function is called by the integrator on state data before using it to evaluate a right-hand side.
    */
    std::function<void (T&, amrex::Real, int, int)> post_update;

    SRI::ButcherTableauTypes butcher_type;
    int number_nodes;
    amrex::Vector<std::unique_ptr<T> > F_nodes;
    amrex::Vector<amrex::Vector<amrex::Real> > tableau;
    amrex::Vector<amrex::Real> weights;
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

    void set_rhs (std::function<void(T&, const T&, T&, const amrex::Real, int)> F)
    {
        rhs = F;
    }

    void set_post_update (std::function<void (T&, amrex::Real, int, int)> F)
    {
        post_update = F;
    }

    std::function<void (T&, amrex::Real, int, int)> get_post_update ()
    {
        return post_update;
    }

    std::function<void(T&, const T&, T&, const amrex::Real, int)> get_rhs ()
    {
        return rhs;
    }

    amrex::Real advance (T& S_old, T& S_new, amrex::Real time, const amrex::Real time_step)
    {
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
        const amrex::Vector<int> scomp_all = {0,0,0,0,0,0,0};
        const amrex::Vector<int> ncomp_all = {nav,1,1,1,nav,nav,nav};

        amrex::Vector<amrex::MultiFab> S_scratch;

        // Fill the RHS F_nodes at each stage
        for (int i = 0; i < number_nodes; ++i)
        {
            // Get current stage time, t = t_old + h * Ci
            amrex::Real stage_time = time + time_step * nodes[i];

            // Fill S_new with the solution value for evaluating F at the current stage
            // Copy S_new = S_old
            amrex::IntegratorOps<T>::Copy(S_new, S_old);
            if (i > 0) {
                // Saxpy across the tableau row:
                // S_new += h * Aij * Fj
                // We should fuse these kernels ...
                for (int j = 0; j < i; ++j)
                {
                    amrex::IntegratorOps<T>::Saxpy(S_new, time_step * tableau[i][j], *F_nodes[j], scomp_all, ncomp_all);
                }

                // Call the post-update hook for the stage state value
                post_update(S_new, stage_time, S_new[IntVar::cons].nGrow(), S_new[IntVar::xmom].nGrow());
            }

            // Fill F[i], the RHS at the current stage
            // F[i] = RHS(y, t) at y = stage_value, t = stage_time
            rhs(*F_nodes[i], S_new, S_scratch, stage_time, RHSVar::all);
        }

        // Fill new State, starting with S_new = S_old.
        // Then Saxpy S_new += h * Wi * Fi for integration weights Wi
        // We should fuse these kernels ...
        amrex::IntegratorOps<T>::Copy(S_new, S_old);
        for (int i = 0; i < number_nodes; ++i)
        {
            amrex::IntegratorOps<T>::Saxpy(S_new, time_step * weights[i], *F_nodes[i], scomp_all, ncomp_all);
        }

        // Call the post-update hook for S_new
        post_update(S_new, time + time_step, S_new[IntVar::cons].nGrow(), S_new[IntVar::xmom].nGrow());

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
