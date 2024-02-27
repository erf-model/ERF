#ifndef MATPROP_NACL_H
#define MATPROP_NACL_H

#include "MaterialProperties.H"

/*! Material properties of NaCl */
class MatProp_NaCl : public MaterialProperties {
    public:

        MatProp_NaCl ()
        {
            m_name = "NaCl";
            m_density = 2170.0; // kg m^{-3}
        }

        virtual ~MatProp_NaCl () = default;

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
