#ifndef _PROBCOMMON_H_
#define _PROBCOMMON_H_

#include <AMReX_ParmParse.H>
#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_MultiFab.H>

#include "DataStruct.H"
#include "EOS.H"
#include "HSE_utils.H"
#include "TileNoZ.H"

struct ProbParmDefaults {
    amrex::Real rho_0 = 1.0;
    amrex::Real T_0 = 300.0;
};

/**
 * Class to hold problem-specific routines
*/
class ProblemBase
{
public:

    /**
     * Virtual destructor to avoid data leakage with derived class
     */
    virtual ~ProblemBase () = default;

    /**
     * Function to initialize the hydrostatic reference density
     *
     * @param[out] rho_hse hydrostatic reference density
     * @param[in] z_phys_nd height coordinate at nodes
     * @param[in] z_phys_cc height coordinate at cell centers
     * @param[in] geom container for geometric information
    */
    virtual void
    erf_init_dens_hse (amrex::MultiFab& /*rho_hse*/,
                       std::unique_ptr<amrex::MultiFab>& /*z_phys_nd*/,
                       std::unique_ptr<amrex::MultiFab>& /*z_phys_cc*/,
                       amrex::Geometry const& /*geom*/)
    {
        amrex::Print() << "Hydrostatically balanced density was NOT set"
            << " -- an appropriate init_type should probably have been specified"
            << " (e.g., input_sounding, ideal, real, or metgrid)"
            << std::endl;
        amrex::Error("Should never call erf_init_dens_hse for "+name()+" problem");
    }

    virtual void
    erf_init_dens_hse_moist (amrex::MultiFab& /*rho_hse*/,
                             std::unique_ptr<amrex::MultiFab>& /*z_phys_nd*/,
                             amrex::Geometry const& /*geom*/)
    {

    }

    /**
     * Function to perform custom initialization of a test problem
     *
     * @param[in]  bx cell-centered box on which to initialize scalars
     * @param[in] xbx box on which to initialize x-component of velocity
     * @param[in] ybx box on which to initialize y-component of velocity
     * @param[in] zbx box on which to initialize z-component of velocity
     * @param[out] state cell-centered variables to be filled in this routine
     * @param[out] x_velocity  x-component of velocity to be filled in this routine
     * @param[out] y_velocity  y-component of velocity to be filled in this routine
     * @param[out] z_velocity  z-component of velocity to be filled in this routine
     * @param[out] r_hse       hydrostatic reference density
     * @param[out] p_hse       hydrostatic reference pressure
     * @param[in] z_nd         height coordinate at nodes
     * @param[in] z_cc         height coordinate at cell centers
     * @param[in] qv           water vapor
     * @param[in] qc           cloud water
     * @param[in] qi           cloud ice
     * @param[in] mf_m         map factor on cell centers
     * @param[in] mf_u         map factor on x-faces
     * @param[in] mf_v         map factor on y-faces
     * @param[in] sc           SolverChoice structure that carries parameters
    */
    virtual void
    init_custom_pert (const amrex::Box&  /*bx*/,
                      const amrex::Box& /*xbx*/,
                      const amrex::Box& /*ybx*/,
                      const amrex::Box& /*zbx*/,
                      amrex::Array4<amrex::Real      > const& /*state*/,
                      amrex::Array4<amrex::Real      > const& /*x_vel*/,
                      amrex::Array4<amrex::Real      > const& /*y_vel*/,
                      amrex::Array4<amrex::Real      > const& /*z_vel*/,
                      amrex::Array4<amrex::Real      > const& /*r_hse*/,
                      amrex::Array4<amrex::Real      > const& /*p_hse*/,
                      amrex::Array4<amrex::Real const> const& /*z_nd*/,
                      amrex::Array4<amrex::Real const> const& /*z_cc*/,
                      amrex::GeometryData const& /*geomdata*/,
                      amrex::Array4<amrex::Real const> const& /*mf_m*/,
                      amrex::Array4<amrex::Real const> const& /*mf_u*/,
                      amrex::Array4<amrex::Real const> const& /*mf_v*/,
                      const SolverChoice& /*sc*/
    )
    {
        amrex::Print() << "No perturbation to background fields supplied for "
            << name() << " problem" << std::endl;;
    }

    /**
     * Function to perform custom initialization of terrain
     *
     * Note: Terrain functionality can also be used to provide grid stretching.
     *
     * @param[in] geom container for geometric information
     * @param[out] z_phys_nd height coordinate at nodes
     * @param[in] time current time
    */
    virtual void
    init_custom_terrain (const amrex::Geometry& /*geom*/,
                         amrex::MultiFab& z_phys_nd,
                         const amrex::Real& /*time*/)
    {
        amrex::Print() << "Initialized flat terrain at z=0" << std::endl;
        z_phys_nd.setVal(0.0);
    }

#ifdef ERF_USE_TERRAIN_VELOCITY
    virtual amrex::Real compute_terrain_velocity(const amrex::Real /*time*/)
    {
        amrex::Error("Should never call compute_terrain_velocity for "+name()+" problem");
    }
#endif

    /**
     * Function to define the quantities needed to impose Rayleigh damping
     *
     * @param[in] tau  strength of Rayleigh damping
     * @param[in] ubar reference value for x-velocity used to define Rayleigh damping
     * @param[in] vbar reference value for y-velocity used to define Rayleigh damping
     * @param[in] wbar reference value for z-velocity used to define Rayleigh damping
     * @param[in] thetabar reference value for potential temperature used to define Rayleigh damping
     * @param[in] geom container for geometric information
    */
    virtual void
    erf_init_rayleigh (amrex::Vector<amrex::Real>& /*tau*/,
                       amrex::Vector<amrex::Real>& /*ubar*/,
                       amrex::Vector<amrex::Real>& /*vbar*/,
                       amrex::Vector<amrex::Real>& /*wbar*/,
                       amrex::Vector<amrex::Real>& /*thetabar*/,
                       amrex::Geometry      const& /*geom*/)
    {
        amrex::Error("Should never call erf_init_rayleigh for "+name()+" problem");
    }

    /**
     * Function to set uniform background density and pressure fields
    */
    void
    init_uniform (const amrex::Box& bx, amrex::Array4<amrex::Real> const& state)
    {
        amrex::Real rho_0 = base_parms.rho_0;
        amrex::Real T_0 = base_parms.T_0;
        amrex::Print() << "Initializing uniform fields"
            << " rho=" << rho_0 << " theta=" << T_0
            << " -- this probably only makes sense with gravity turned off"
            << std::endl;

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            state(i, j, k, Rho_comp) = rho_0;
            state(i, j, k, RhoTheta_comp) = rho_0 * T_0;
        });
    }

protected:
    // Struct to store problem parameters
    ProbParmDefaults base_parms;

    /**
     * Function to update default base parameters, currently only used for
     * init_type=='uniform'
     */
    void init_base_parms (amrex::Real rho_0, amrex::Real T_0) {
        base_parms.rho_0 = rho_0;
        base_parms.T_0 = T_0;
    }

    // descriptor for problem definition
    virtual std::string name () = 0;
};


/**
 * Function to init the physical bounds of the domain
 * and instantiate a Problem derived from ProblemBase
*/
std::unique_ptr<ProblemBase> amrex_probinit (const amrex_real* problo,
                                             const amrex_real* probhi);

#endif
