#ifndef MATERIALPROPERTIES_H
#define MATERIALPROPERTIES_H

#include <string>
#include <AMReX_REAL.H>
#include <AMReX_Gpu.H>

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

/*! Class for material properties */
class MaterialProperties {
    public:

        MaterialProperties ( const std::string& a_name );

        ~MaterialProperties () = default;

        inline amrex::Real density() const { return m_density; }
        inline const std::string& name() const { return m_name; }
        inline amrex::Real term_vel_a() const { return m_a_tv; }
        inline amrex::Real term_vel_b() const { return m_b_tv; }

    protected:

        std::string m_name;       /*!< name */
        amrex::Real m_density;    /*!< density */

        /* terminal velocity coeffs */
        amrex::Real m_a_tv, m_b_tv;

    private:

        void setProperties_H2O();
        void setProperties_NaCl();

};

#endif
