#ifndef MATERIALPROPERTIES_H
#define MATERIALPROPERTIES_H

#include <AMReX_Enum.H>
#include <AMReX_REAL.H>
#include <AMReX_Gpu.H>
#include <AMReX_MultiFab.H>
#include "ERF_Constants.H"

/*! Materials */
namespace Species
{
    AMREX_ENUM(Name,
        H2O,
        NaCl,
        NH42SO4,
        NH4HSO4,
        soil,
        water, // same as H2O, to test multispecies SDM
        agua,  // same as H2O, to test multispecies SDM
        none
    );
}

/*! Functions for computing saturation pressure and vapour fraction */
namespace saturation_funcs
{
    /*! \brief Null function for saturation pressure */
    AMREX_GPU_HOST
    void compute_saturation_pressure_null (  amrex::MultiFab&, const amrex::MultiFab& );
    /*! \brief Compute saturation pressure for moist air */
    AMREX_GPU_HOST
    void compute_saturation_pressure_H2O  (  amrex::MultiFab&, const amrex::MultiFab& );

    /*! \brief Null function for saturation vapour fraction */
    AMREX_GPU_HOST
    void compute_saturation_vapfrac_null ( amrex::MultiFab&,
                                           const amrex::MultiFab&,
                                           const amrex::MultiFab& );
    /*! \brief Compute saturation vapour fraction for moist air */
    AMREX_GPU_HOST
    void compute_saturation_vapfrac_H2O  ( amrex::MultiFab&,
                                           const amrex::MultiFab&,
                                           const amrex::MultiFab& );
}

/*! Class for material properties */
struct MaterialProperties {

    static constexpr amrex::Real s_N_av = 6.02214076e23; /*!< Avogadro's number */
    static constexpr amrex::Real s_kb   = 1.380649e-23; /*!< Boltzmann constant [J K^{-1}] */

    /*! Molecular weight of air [kg mol^{-1}] */
    static constexpr amrex::Real s_mol_weight_air = 0.0289647;
    /*! Molecular length scale of air [m] */
    static constexpr amrex::Real s_sigma_air = 3.62e-10;
    /*! Temperature parameter for air [K] */
    static constexpr amrex::Real s_eps_air_k = 97;

    Species::Name m_name = Species::Name::none; /*!< name */
    amrex::Real m_density = DBL_MAX; /*!< density */
    amrex::Real m_ionization = DBL_MAX; /*!< ionization */
    amrex::Real m_mol_weight = DBL_MAX; /*!< molecular weight (condensate) */
    amrex::Real m_lat_vap = DBL_MAX; /*!< latent heat of vaporization */
    amrex::Real m_Rv = DBL_MAX; /*!< gas constant for vapour of this material */
    amrex::Real m_Tc = DBL_MAX; /*!< critical temperature for surface tension */
    amrex::Real m_Tb = DBL_MAX; /*!< boiling temperature */
    amrex::Real m_Nav_by_molweight = DBL_MAX; /*!< Avogadro number by vapour mol. weight */
    /*! vector of molar heat capacity polynomial coeffs */
    amrex::Real m_mol_Cp_coeffs[7] = {  amrex::Real(DBL_MAX),
                                        amrex::Real(DBL_MAX),
                                        amrex::Real(DBL_MAX),
                                        amrex::Real(DBL_MAX),
                                        amrex::Real(DBL_MAX),
                                        amrex::Real(DBL_MAX),
                                        amrex::Real(DBL_MAX) };
    bool m_is_soluble = false; /*!< is soluble in water? */
    bool m_is_water = false; /*!< is this water? */
    /*! pointer to function to compute saturation pressure */
    decltype(saturation_funcs::compute_saturation_pressure_null)* m_saturation_pressure_func = nullptr;
    /*! pointer to function to compute vapour fraction */
    decltype(saturation_funcs::compute_saturation_vapfrac_null)* m_saturation_vapfrac_func = nullptr;

    /*! \brief Constructor */
    AMREX_GPU_HOST_DEVICE
    MaterialProperties ( const Species::Name& a_name );

    /*! \brief Copy constructor */
    AMREX_GPU_HOST_DEVICE
    MaterialProperties ( const MaterialProperties& a_matprop );

    /*! \brief Default destructor */
    AMREX_GPU_HOST_DEVICE
    ~MaterialProperties () = default;

    /*! \brief Print parameters */
    AMREX_GPU_HOST
    void print() const {
        amrex::Print() << "Material properties of " << amrex::getEnumNameString(m_name) << ":\n";
        amrex::Print() << "    density: " << m_density << "\n";
        amrex::Print() << "    ionization: " << m_ionization << "\n";
        amrex::Print() << "    mol. weight: " << m_mol_weight << "\n";
        amrex::Print() << "    latent heat (vap.): " << m_lat_vap << "\n";
        amrex::Print() << "    Rv: " << m_Rv << "\n";
        amrex::Print() << "    Tc: " << m_Tc << "\n";
        amrex::Print() << "    Tb: " << m_Tb << "\n";
        amrex::Print() << "    N_av/mol.weight: " << m_Nav_by_molweight << "\n";
        amrex::Print() << "    mol. Cp coeffs: ";
        for (int i = 0; i < 7; i++) { amrex::Print() << m_mol_Cp_coeffs[i] << ", "; }
        amrex::Print() << "\n";
    }

    /*! \brief Set this material to H20 */
    AMREX_GPU_HOST_DEVICE
    void setProperties_H2O();
    /*! \brief Set this material to water */
    AMREX_GPU_HOST_DEVICE
    void setProperties_water();
    /*! \brief Set this material to agua */
    AMREX_GPU_HOST_DEVICE
    void setProperties_agua();
    /*! \brief Set this material to NaCl */
    AMREX_GPU_HOST_DEVICE
    void setProperties_NaCl();
    /*! \brief Set this material to ammonium sulfate */
    AMREX_GPU_HOST_DEVICE
    void setProperties_NH42SO4();
    /*! \brief Set this material to ammonium bisulfate */
    AMREX_GPU_HOST_DEVICE
    void setProperties_NH4HSO4();
    /*! \brief Set this material to soil */
    AMREX_GPU_HOST_DEVICE
    void setProperties_soil();

    /*! \brief Return the molar heat capacity */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real molarHeatCapacity( const amrex::Real a_T ) const
    {
        auto retval = 0.0;
        auto term = 1.0;
        for (auto i = 0; i < 7; i++) {
            retval += (m_mol_Cp_coeffs[i] * term);
            term *= a_T;
        }
        return retval;
    }

    /*! \brief Compute the surface tension coeff given temperature */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real surfaceTension( const amrex::Real a_T ) const
    {
        if (m_is_water) {
            return amrex::Real(0.076148325); // N m^{-1}
        } else {
            auto term1 = 6.0 * std::cbrt(s_N_av);
            auto term2 = std::cbrt((m_density/m_mol_weight)*(m_density/m_mol_weight));
            auto retval = molarHeatCapacity(a_T)*(term2/term1)*(m_Tc - a_T);
            return retval;
        }
    }

    /*! \brief Return the coeff of curvature given the temperature */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real coeffCurv( const amrex::Real a_T ) const
    {
        auto surf_tens = surfaceTension(a_T);
        auto retval = 2 * surf_tens / (m_Rv * m_density);
        return retval;
    }

    /*! \brief Return the vapour pressure coefficient */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real coeffVPSolute() const
    {
        auto retval = (0.75/PI)*(m_mol_weight/m_density);
        if (m_is_water) {
            retval = 4.3e-06;
        }
        return retval;
    }

    /*! \brief Return the molecular diffusion coefficient given temperature and pressure */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real coeffMolecularDiffusion( const amrex::Real a_T, const amrex::Real a_P ) const
    {
        if (m_is_water) {

            return diffelq; // from ERF_Constants.H (2.21e-05 [m^2 s^{-1}])

        } else {

            auto mol_weight_v = s_N_av/m_Nav_by_molweight;
            auto sigma = 1.18e-8 * std::cbrt(mol_weight_v/m_density);
            auto sigma_v_air = 0.5 * (sigma + s_sigma_air);

            auto eps_v = 1.15 * m_Tb * s_kb;
            auto eps_v_air = std::sqrt( eps_v * s_eps_air_k * s_kb);
            auto T_star = a_T / (eps_v_air/s_kb);

            auto omegaD =    1.06036 / std::pow(T_star, 0.1561)
                          +  0.193   / std::exp(0.47635*T_star)
                          +  1.03587 / std::exp(1.52996*T_star)
                          +  1.76474 / std::exp(3.89411*T_star);

            auto term1 = 3 * std::sqrt( s_kb*a_T * s_kb*a_T * s_kb*a_T);
            auto term2 = 8.0 * std::sqrt(2*PI) * a_P * sigma_v_air*sigma_v_air * omegaD;
            auto term3 = std::sqrt( m_Nav_by_molweight + s_N_av/s_mol_weight_air );

            auto retval = term1*term3/term2;
            return retval;
        }
    }

    /*! \brief Compute saturation pressure */
    AMREX_GPU_HOST AMREX_FORCE_INLINE
    void computeSaturationPressure ( amrex::MultiFab& a_e, const amrex::MultiFab& a_T) const
    {
        m_saturation_pressure_func(a_e, a_T);
    }

    /*! \brief Compute saturation vapour fraction */
    AMREX_GPU_HOST AMREX_FORCE_INLINE
    void computeSaturationVapFrac  ( amrex::MultiFab&         a_S,
                                     const amrex::MultiFab&   a_T,
                                     const amrex::MultiFab&   a_p ) const
    {
        m_saturation_vapfrac_func(a_S, a_T, a_p);
    }

};

#endif
