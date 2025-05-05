#include "ERF_MaterialProperties.H"
#include "ERF_MicrophysicsUtils.H"

using namespace amrex;

namespace saturation_funcs
{
    void compute_saturation_pressure_null ( MultiFab&, const MultiFab&) { }

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

    void compute_saturation_vapfrac_null ( MultiFab&, const MultiFab&) { }

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
                                            pressure_arr(i,j,k,0)/100.0,
                                            qsat_arr(i,j,k,0) );
                            } );
        }
    }
}

MaterialProperties::MaterialProperties ( const std::string& a_name )
{
    m_name = a_name;
    if (a_name == MaterialNames::h2o) {
        setProperties_H2O();
    } else if (a_name == MaterialNames::nacl) {
        setProperties_NaCl();
    } else if (a_name == MaterialNames::amsu) {
        setProperties_AmSu();
    } else if (a_name == MaterialNames::soil) {
        setProperties_soil();
    } else {
        amrex::Abort("ERROR: undefined material in MaterialProperties()");
    }
}

void MaterialProperties::setProperties_H2O()
{
    m_density = rhor; // ERF_Constants.H

    m_coeff_curv = 3.3e-07; // m K
    m_coeff_VP_solute = 4.3e-06; // m^3
    m_ionization = 2;
    m_mol_weight = 1.802e-02; // kg mol^-1
    m_lat_vap = L_v; // ERF_Constants.H
    m_Rv = R_v; // ERF_Constants.H
    m_is_soluble = true;

    m_saturation_pressure_func = saturation_funcs::compute_saturation_pressure_H2O;
    m_saturation_vapfrac_func = saturation_funcs::compute_saturation_vapfrac_H2O;

}

void MaterialProperties::setProperties_NaCl()
{
    m_density = 2170.0;

    m_coeff_curv = DBL_MAX; // m K
    m_coeff_VP_solute = DBL_MAX; // m^3
    m_ionization = 2;
    m_mol_weight = 5.844e-02; //kg mol^-1
    m_lat_vap = DBL_MAX;
    m_Rv = DBL_MAX;
    m_is_soluble = true;

    m_saturation_pressure_func = nullptr;
    m_saturation_vapfrac_func = nullptr;
}

void MaterialProperties::setProperties_AmSu()
{
    m_density = 1780.0;

    m_coeff_curv = DBL_MAX; // m K
    m_coeff_VP_solute = DBL_MAX; // m^3
    m_ionization = 2;
    m_mol_weight = 1.11511e-01; //kg mol^-1
    m_lat_vap = DBL_MAX;
    m_Rv = DBL_MAX;
    m_is_soluble = true;

    m_saturation_pressure_func = nullptr;
    m_saturation_vapfrac_func = nullptr;
}

void MaterialProperties::setProperties_soil()
{
    m_density = 1140.0; // NNSS Area 5 sample (Spriggs and Ray-Maitra, 2007)

    m_coeff_curv = DBL_MAX; // m K
    m_coeff_VP_solute = DBL_MAX; // m^3
    m_ionization = DBL_MAX;
    m_mol_weight = DBL_MAX; //kg mol^-1
    m_lat_vap = DBL_MAX;
    m_Rv = DBL_MAX;
    m_is_soluble = false;

    m_saturation_pressure_func = nullptr;
    m_saturation_vapfrac_func = nullptr;
}
