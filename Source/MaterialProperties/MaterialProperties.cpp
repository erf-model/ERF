#include "MaterialProperties.H"

MaterialProperties::MaterialProperties ( const std::string& a_name )
{
    m_name = a_name;
    if (a_name == MaterialNames::h2o) {
        setProperties_H2O();
    } else if (a_name == MaterialNames::nacl) {
        setProperties_NaCl();
    } else {
        amrex::Abort("ERROR: undefined material in MaterialProperties()");
    }
}

void MaterialProperties::setProperties_H2O()
{
    m_density = 1.0;

    m_a_tv = 3.778;
    m_b_tv = 0.67;
}

void MaterialProperties::setProperties_NaCl()
{
    m_density = 2170.0;

    m_a_tv = 3.778;
    m_b_tv = 0.67;
}
