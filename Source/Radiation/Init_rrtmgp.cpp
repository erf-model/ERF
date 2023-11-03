#include "mo_gas_concentrations.h"
#include "mo_gas_optics_rrtmgp.h"
#include "mo_load_coefficients.h"
#include "mo_rte_sw.h"
#include "mo_rte_lw.h"
#include "mo_optical_props.h"
#include "rrtmgp_const.h"
#include "mo_fluxes_byband.h"

// Rrtmgp
#include "Rrtmgp.H"

void Rrtmgp::initialize(int num_gas, std::vector<std::string> active_gas_names,
                        const char* rrtmgp_coefficients_file_sw,
                        const char* rrtmgp_coefficients_file_lw)
{
    // First, make sure yakl has been initialized
    if (!yakl::isInitialized()) {
        yakl::init();
    }

    // Read gas optics coefficients from file
    // Need to initialize available_gases here! The only field of the
    // available_gases type that is used in the kdist initialize is
    // available_gases%gas_name, which gives the name of each gas that would be
    // present in the ty_gas_concs object. So, we can just set this here, rather
    // than trying to fully populate the ty_gas_concs object here, which would be
    // impossible from this initialization routine because I do not think the
    // rad_cnst objects are setup yet.
    // the other tasks!
    ngas      = num_gas;
    for (auto i = 0; i < ngas; ++i)
       gas_names[i] = active_gas_names[i].c_str();

    coefficients_file_sw = rrtmgp_coefficients_file_sw;
    coefficients_file_lw = rrtmgp_coefficients_file_lw;

    auto active_gases = string1d("active_gases", ngas);
    for (int igas=0; igas<ngas; igas++) {
        active_gases(igas+1) = gas_names[igas];
    }
    GasConcs available_gases;
    available_gases.init(active_gases, 1, 1);
    load_and_init(k_dist_sw, coefficients_file_sw, available_gases);
    load_and_init(k_dist_lw, coefficients_file_lw, available_gases);
}

