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

void Rrtmgp::finalize ()
{
    k_dist_sw.finalize();
    k_dist_lw.finalize();
    yakl::finalize();
}

