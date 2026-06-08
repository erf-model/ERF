#include "ERF_MaterialProperties.H"
#include "ERF_MicrophysicsUtils.H"

using namespace amrex;

namespace saturation_funcs
{
    AMREX_GPU_HOST
    void compute_saturation_pressure_null ( MultiFab&, const MultiFab&) { }

    AMREX_GPU_HOST
    void compute_saturation_pressure_H2O  ( MultiFab&       a_mf_sat_pressure,
                                            const MultiFab& a_mf_temperature)
    {
        const auto& gvec = a_mf_sat_pressure.nGrowVect();
        for (MFIter mfi(a_mf_sat_pressure, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            Box bx = mfi.tilebox();
            bx.grow(gvec);
            const Array4<Real>& psat_arr = a_mf_sat_pressure.array(mfi);
            const Array4<Real const>& temperature_arr = a_mf_temperature.array(mfi);

            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                            { psat_arr(i,j,k,0) = erf_esatw(temperature_arr(i,j,k,0))*100; } );
                              // formula gives pressure in hPa; we will save it in Pa.
        }
    }

    AMREX_GPU_HOST
    void compute_saturation_vapfrac_null ( MultiFab&, const MultiFab&) { }

    AMREX_GPU_HOST
    void compute_saturation_vapfrac_H2O ( MultiFab&          a_mf_sat_vapfrac,
                                          const MultiFab&    a_mf_temperature,
                                          const MultiFab&    a_mf_pressure )
    {
        const auto& gvec = a_mf_sat_vapfrac.nGrowVect();
        for (MFIter mfi(a_mf_sat_vapfrac, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            Box bx = mfi.tilebox();
            bx.grow(gvec);
            const Array4<Real>& qsat_arr = a_mf_sat_vapfrac.array(mfi);
            const Array4<Real const>& temperature_arr = a_mf_temperature.array(mfi);
            const Array4<Real const>& pressure_arr = a_mf_pressure.array(mfi);

            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                            {
                                // pressure is in Pa; formula takes pressure in hPa
                                erf_qsatw(  temperature_arr(i,j,k,0),
                                            pressure_arr(i,j,k,0)/Real(100.0),
                                            qsat_arr(i,j,k,0) );
                            } );
        }
    }
}

AMREX_GPU_HOST_DEVICE
MaterialProperties::MaterialProperties ( const Species::Name& a_name )
{
    m_name = a_name;

    if (a_name == Species::Name::H2O) {
        setProperties_H2O();
    } else if (a_name == Species::Name::water) {
        setProperties_water();
    } else if (a_name == Species::Name::agua) {
        setProperties_agua();
    } else if (a_name == Species::Name::NaCl) {
        setProperties_NaCl();
    } else if (a_name == Species::Name::NH42SO4) {
        setProperties_NH42SO4();
    } else if (a_name == Species::Name::NH4HSO4) {
        setProperties_NH4HSO4();
    } else if (a_name == Species::Name::soil) {
        setProperties_soil();
    } else {
        amrex::Abort("ERROR: undefined material in MaterialProperties()");
    }
    m_is_soluble = (m_ionization > 0);
}

AMREX_GPU_HOST_DEVICE
MaterialProperties::MaterialProperties ( const MaterialProperties& a_matprop )
{
    m_name = a_matprop.m_name;
    m_density = a_matprop.m_density;
    m_ionization = a_matprop.m_ionization;
    m_mol_weight = a_matprop.m_mol_weight;
    m_lat_vap = a_matprop.m_lat_vap;
    m_Rv = a_matprop.m_Rv;
    m_Tc = a_matprop.m_Tc;
    m_Tb = a_matprop.m_Tb;
    m_Nav_by_molweight = a_matprop.m_Nav_by_molweight;
    for (auto i = 0; i < 7; i++) { m_mol_Cp_coeffs[i] = a_matprop.m_mol_Cp_coeffs[i]; }
    m_is_soluble = a_matprop.m_is_soluble;
    m_is_water = a_matprop.m_is_water;
    AMREX_IF_ON_HOST((
        m_saturation_pressure_func = a_matprop.m_saturation_pressure_func;
        m_saturation_vapfrac_func = a_matprop.m_saturation_vapfrac_func;
    ))
}

AMREX_GPU_HOST_DEVICE
void MaterialProperties::setProperties_H2O()
{
    m_density = rhor; // ERF_Constants.H

    m_ionization = 0;
    m_mol_weight = Real(1.802e-02); // kg mol^-1
    m_lat_vap = L_v; // ERF_Constants.H
    m_Rv = R_v; // ERF_Constants.H
    m_Tb = Real(373.15); // K
    m_Nav_by_molweight = s_N_av / m_mol_weight;
    m_is_water = true;

    AMREX_IF_ON_HOST((
        m_saturation_pressure_func = saturation_funcs::compute_saturation_pressure_H2O;
        m_saturation_vapfrac_func = saturation_funcs::compute_saturation_vapfrac_H2O;
    ))

}

AMREX_GPU_HOST_DEVICE
void MaterialProperties::setProperties_water()
{
    setProperties_H2O();
}

AMREX_GPU_HOST_DEVICE
void MaterialProperties::setProperties_agua()
{
    setProperties_H2O();
}

AMREX_GPU_HOST_DEVICE
void MaterialProperties::setProperties_NaCl()
{
    m_density = Real(2170.0);

    m_ionization = 2;
    m_mol_weight = Real(5.844e-02); //kg mol^-1

    m_saturation_pressure_func = nullptr;
    m_saturation_vapfrac_func = nullptr;
}

AMREX_GPU_HOST_DEVICE
void MaterialProperties::setProperties_NH42SO4()
{
    m_density = Real(1770.0);

    m_ionization = 3; // 2xNH4 + 1xSO4
    m_mol_weight = Real(1.3214e-01); //kg mol^-1

    m_saturation_pressure_func = nullptr;
    m_saturation_vapfrac_func = nullptr;
}

AMREX_GPU_HOST_DEVICE
void MaterialProperties::setProperties_NH4HSO4()
{
    m_density = Real(1780.0);

    m_ionization = 2; // NH4+ and HSO4-
    m_mol_weight = Real(1.1511e-01); //kg mol^-1

    m_saturation_pressure_func = nullptr;
    m_saturation_vapfrac_func = nullptr;
}

AMREX_GPU_HOST_DEVICE
void MaterialProperties::setProperties_soil()
{
    m_density = Real(1220.0); // loose dry dirt

    m_ionization = 0;

    m_saturation_pressure_func = nullptr;
    m_saturation_vapfrac_func = nullptr;
}
