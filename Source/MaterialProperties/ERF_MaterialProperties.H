#ifndef MATERIALPROPERTIES_H
#define MATERIALPROPERTIES_H

#include <string>
#include <AMReX_REAL.H>
#include <AMReX_Gpu.H>
#include <AMReX_MultiFab.H>
#include "ERF_Constants.H"

/*! Material names */
namespace MaterialNames
{
    /*! Water */
    const std::string h2o = "H2O";
    /*! also water */
    const std::string water = "water";
    /*! also water */
    const std::string agua  = "agua";
    /*! Salt (NaCl) */
    const std::string nacl = "NaCl";
    /*! Ammonium Sulfate/Bisulfate */
    const std::string amsu = "NH42SO4";
    /*! Soil */
    const std::string soil = "soil";
}

/*! Functions for computing saturation pressure and vapour fraction */
namespace saturation_funcs
{
    /*! \brief Null function for saturation pressure */
    void compute_saturation_pressure_null (  amrex::MultiFab&, const amrex::MultiFab& );
    /*! \brief Compute saturation pressure for moist air */
    void compute_saturation_pressure_H2O  (  amrex::MultiFab&, const amrex::MultiFab& );

    /*! \brief Null function for saturation vapour fraction */
    void compute_saturation_vapfrac_null ( amrex::MultiFab&,
                                           const amrex::MultiFab&,
                                           const amrex::MultiFab& );
    /*! \brief Compute saturation vapour fraction for moist air */
    void compute_saturation_vapfrac_H2O  ( amrex::MultiFab&,
                                           const amrex::MultiFab&,
                                           const amrex::MultiFab& );
}

/*! Class for material properties */
class MaterialProperties {
    public:

        /*! \brief Constructor */
        MaterialProperties ( const std::string& a_name );

        /*! \brief Default destructor */
        ~MaterialProperties () = default;

        /*! \brief Return the material density */
        [[nodiscard]] inline amrex::Real density() const { return m_density; }
        /*! \brief Return the material name */
        [[nodiscard]] inline const std::string& name() const { return m_name; }

        /*! \brief Return the coeff of curvature */
        [[nodiscard]] inline amrex::Real coeffCurv() const { return m_coeff_curv; }
        /*! \brief Return the vapour pressure coefficient */
        [[nodiscard]] inline amrex::Real coeffVPSolute() const { return m_coeff_VP_solute*m_ionization; }
        /*! \brief Return the molecular weight */
        [[nodiscard]] inline amrex::Real molWeight() const { return m_mol_weight; }
        /*! \brief Return the latent heat of vaporization */
        [[nodiscard]] inline amrex::Real latHeatVap() const { return m_lat_vap; }
        /*! \brief Return the thermal conductivity */
        [[nodiscard]] inline amrex::Real thermCond() const { return m_therm_cond; }
        /*! \brief Return the gas constant for air with vapour of this material */
        [[nodiscard]] inline amrex::Real Rv() const { return m_Rv; }

        /*! \brief Returns whether material is soluble */
        [[nodiscard]] inline bool isSoluble() const { return m_is_soluble; }

        /*! \brief Compute saturation pressure */
        inline void computeSaturationPressure ( amrex::MultiFab& a_e, const amrex::MultiFab& a_T) const
        {
            m_saturation_pressure_func(a_e, a_T);
        }

        /*! \brief Compute saturation vapour fraction */
        inline void computeSaturationVapFrac  ( amrex::MultiFab&         a_S,
                                                const amrex::MultiFab&   a_T,
                                                const amrex::MultiFab&   a_p ) const
        {
            m_saturation_vapfrac_func(a_S, a_T, a_p);
        }

    protected:

        std::string m_name; /*!< name */
        amrex::Real m_density; /*!< density */

        /* various coeffs */
        amrex::Real m_coeff_curv; /*!< curvature coefficient for condensation/evaporation equation */
        amrex::Real m_coeff_VP_solute; /*!< vapor pressure coefficient for condensation/evaporation equation */
        amrex::Real m_ionization; /*!< ionization */
        amrex::Real m_mol_weight; /*!< molecular weight */
        amrex::Real m_lat_vap; /*!< latent heat of vaporization */
        amrex::Real m_therm_cond; /*!< thermal conductivity */
        amrex::Real m_Rv; /*!< gas constant for vapour of this material */

        bool m_is_soluble; /*!< is soluble in water? */

        /*! pointer to function to compute saturation pressure */
        decltype(saturation_funcs::compute_saturation_pressure_null)* m_saturation_pressure_func;
        /*! pointer to function to compute vapour fraction */
        decltype(saturation_funcs::compute_saturation_vapfrac_null)* m_saturation_vapfrac_func;

    private:

        /*! \brief Set this material to H20 */
        void setProperties_H2O();
        /*! \brief Set this material to water */
        void setProperties_water();
        /*! \brief Set this material to agua */
        void setProperties_agua();
        /*! \brief Set this material to NaCl */
        void setProperties_NaCl();
        /*! \brief Set this material to ammonium sulfate */
        void setProperties_AmSu();
        /*! \brief Set this material to soil */
        void setProperties_soil();

};

#endif
