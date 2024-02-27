#ifndef MATPROP_WATER_H
#define MATPROP_WATER_H

#include "MaterialProperties.H"

/*! Material properties of water */
class MatProp_Water : public MaterialProperties {
    public:

        MatProp_Water ()
        {
            m_name = "H2O";
            m_density = 1.0; // kg m^{-3}
        }

        virtual ~MatProp_Water () = default;

        virtual inline amrex::Real density() const override
        {
            return m_density;
        }

        virtual inline const std::string& name() const
        {
            return m_name;
        }

    protected:

        std::string m_name;       /*!< name */
        amrex::Real m_density;    /*!< density */

    private:

};

#endif
