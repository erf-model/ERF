#ifndef MATERIALPROPERTIES_H
#define MATERIALPROPERTIES_H

#include <string>
#include <AMReX_REAL.H>
#include <AMReX_Gpu.H>
#include <AMReX_MultiFab.H>

namespace MaterialNames
{
    const std::string h2o = "H2O";
    const std::string nacl = "NaCl";
}

namespace MatPropFunctions
{
    struct TerminalVelocity
    {
        amrex::Real a;
        amrex::Real b;

        AMREX_GPU_HOST_DEVICE
        amrex::Real operator() (const amrex::Real a_r) const noexcept
        {
            amrex::Real D = 2*a_r*1000; // meter -> mm
            return a*(std::exp(b*std::log(D)));
        }
    };
}

namespace saturation_funcs
{
    void compute_saturation_pressure_null (  amrex::MultiFab&, const amrex::MultiFab& );
    void compute_saturation_pressure_H2O  (  amrex::MultiFab&, const amrex::MultiFab& );

    void compute_saturation_vapfrac_null ( amrex::MultiFab&,
                                           const amrex::MultiFab&,
                                           const amrex::MultiFab& );
    void compute_saturation_vapfrac_H2O  ( amrex::MultiFab&,
                                           const amrex::MultiFab&,
                                           const amrex::MultiFab& );
}

/*! Class for material properties */
class MaterialProperties {
    public:

        MaterialProperties ( const std::string& a_name );

        ~MaterialProperties () = default;

        inline amrex::Real density() const { return m_density; }
        inline const std::string& name() const { return m_name; }
        inline amrex::Real term_vel_a() const { return m_a_tv; }
        inline amrex::Real term_vel_b() const { return m_b_tv; }

        inline amrex::Real coeffCurv() const { return m_coeff_curv; }
        inline amrex::Real coeffVPSolute(const MaterialProperties& a_solute) const
        {
            return m_coeff_VP_solute * m_ionization / a_solute.molWeight();
        }
        inline amrex::Real molWeight() const { return m_mol_weight; }
        inline amrex::Real latHeatVap() const { return m_lat_vap; }
        inline amrex::Real thermCond() const { return m_therm_cond; }
        inline amrex::Real Rv() const { return m_Rv; }
        inline amrex::Real molDiffCoeff() const { return m_mol_diff; }

        inline void computeSaturationPressure ( amrex::MultiFab& a_e, const amrex::MultiFab& a_T)
        {
            m_saturation_pressure_func(a_e, a_T);
        }
        inline void computeSaturationVapFrac  ( amrex::MultiFab&         a_S,
                                                const amrex::MultiFab&   a_T,
                                                const amrex::MultiFab&   a_p )
        {
            m_saturation_vapfrac_func(a_S, a_T, a_p);
        }

    protected:

        std::string m_name;       /*!< name */
        amrex::Real m_density;    /*!< density */

        /* terminal velocity coeffs */
        amrex::Real m_a_tv, m_b_tv;

        /* other coeffs */
        amrex::Real m_coeff_curv;
        amrex::Real m_coeff_VP_solute;
        amrex::Real m_ionization;
        amrex::Real m_mol_weight;
        amrex::Real m_lat_vap;
        amrex::Real m_therm_cond;
        amrex::Real m_Rv;
        amrex::Real m_mol_diff;

        decltype(saturation_funcs::compute_saturation_pressure_null)* m_saturation_pressure_func;
        decltype(saturation_funcs::compute_saturation_vapfrac_null)* m_saturation_vapfrac_func;

    private:

        void setProperties_H2O();
        void setProperties_NaCl();

};

#endif
